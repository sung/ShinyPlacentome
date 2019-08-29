(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bam.js: indexed binary alignments
//

"use strict";

if (typeof require !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var bin = require('./bin');
    var readInt = bin.readInt;
    var readShort = bin.readShort;
    var readByte = bin.readByte;
    var readInt64 = bin.readInt64;
    var readFloat = bin.readFloat;

    var lh3utils = require('./lh3utils');
    var readVob = lh3utils.readVob;
    var unbgzf = lh3utils.unbgzf;
    var reg2bins = lh3utils.reg2bins;
    var Chunk = lh3utils.Chunk;
}

var BAM_MAGIC = 0x14d4142;
var BAI_MAGIC = 0x1494142;

var BamFlags = {
    MULTIPLE_SEGMENTS: 0x1,
    ALL_SEGMENTS_ALIGN: 0x2,
    SEGMENT_UNMAPPED: 0x4,
    NEXT_SEGMENT_UNMAPPED: 0x8,
    REVERSE_COMPLEMENT: 0x10,
    NEXT_REVERSE_COMPLEMENT: 0x20,
    FIRST_SEGMENT: 0x40,
    LAST_SEGMENT: 0x80,
    SECONDARY_ALIGNMENT: 0x100,
    QC_FAIL: 0x200,
    DUPLICATE: 0x400,
    SUPPLEMENTARY: 0x800
};

function BamFile() {}

// Calculate the length (in bytes) of the BAI ref starting at offset.
// Returns {nbin, length, minBlockIndex}
function _getBaiRefLength(uncba, offset) {
    var p = offset;
    var nbin = readInt(uncba, p);p += 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(uncba, p);
        var nchnk = readInt(uncba, p + 4);
        p += 8 + nchnk * 16;
    }
    var nintv = readInt(uncba, p);p += 4;

    var minBlockIndex = 1000000000;
    var q = p;
    for (var i = 0; i < nintv; ++i) {
        var v = readVob(uncba, q);q += 8;
        if (v) {
            var bi = v.block;
            if (v.offset > 0) bi += 65536;

            if (bi < minBlockIndex) minBlockIndex = bi;
            break;
        }
    }
    p += nintv * 8;

    return {
        minBlockIndex: minBlockIndex,
        nbin: nbin,
        length: p - offset
    };
}

function makeBam(data, bai, indexChunks, callback, attempted) {
    // Do an initial probe on the BAM file to catch any mixed-content errors.
    data.slice(0, 10).fetch(function (header) {
        if (header) {
            return makeBam2(data, bai, indexChunks, callback, attempted);
        } else {
            return callback(null, "Couldn't access BAM.");
        }
    }, { timeout: 5000 });
}

function makeBam2(data, bai, indexChunks, callback, attempted) {
    var bam = new BamFile();
    bam.data = data;
    bam.bai = bai;
    bam.indexChunks = indexChunks;

    var minBlockIndex = bam.indexChunks ? bam.indexChunks.minBlockIndex : 1000000000;

    // Fills out bam.chrToIndex and bam.indexToChr based on the first few bytes of the BAM.
    function parseBamHeader(r) {
        if (!r) {
            return callback(null, "Couldn't access BAM");
        }

        var unc = unbgzf(r, r.byteLength);
        var uncba = new Uint8Array(unc);

        var magic = readInt(uncba, 0);
        if (magic != BAM_MAGIC) {
            return callback(null, "Not a BAM file, magic=0x" + magic.toString(16));
        }
        var headLen = readInt(uncba, 4);
        var header = '';
        for (var i = 0; i < headLen; ++i) {
            header += String.fromCharCode(uncba[i + 8]);
        }

        var nRef = readInt(uncba, headLen + 8);
        var p = headLen + 12;

        bam.chrToIndex = {};
        bam.indexToChr = [];
        for (var i = 0; i < nRef; ++i) {
            var lName = readInt(uncba, p);
            var name = '';
            for (var j = 0; j < lName - 1; ++j) {
                name += String.fromCharCode(uncba[p + 4 + j]);
            }
            var lRef = readInt(uncba, p + lName + 4);
            bam.chrToIndex[name] = i;
            if (name.indexOf('chr') == 0) {
                bam.chrToIndex[name.substring(3)] = i;
            } else {
                bam.chrToIndex['chr' + name] = i;
            }
            bam.indexToChr.push(name);

            p = p + 8 + lName;
        }

        if (bam.indices) {
            return callback(bam);
        }
    }

    function parseBai(header) {
        if (!header) {
            return "Couldn't access BAI";
        }

        var uncba = new Uint8Array(header);
        var baiMagic = readInt(uncba, 0);
        if (baiMagic != BAI_MAGIC) {
            return callback(null, 'Not a BAI file, magic=0x' + baiMagic.toString(16));
        }

        var nref = readInt(uncba, 4);

        bam.indices = [];

        var p = 8;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var o = _getBaiRefLength(uncba, blockStart);
            p += o.length;

            minBlockIndex = Math.min(o.minBlockIndex, minBlockIndex);

            var nbin = o.nbin;

            if (nbin > 0) {
                bam.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
            }
        }

        return true;
    }

    if (!bam.indexChunks) {
        bam.bai.fetch(function (header) {
            // Do we really need to fetch the whole thing? :-(
            var result = parseBai(header);
            if (result !== true) {
                if (bam.bai.url && typeof attempted === "undefined") {
                    // Already attempted x.bam.bai not there so now trying x.bai
                    bam.bai.url = bam.data.url.replace(new RegExp('.bam$'), '.bai');

                    // True lets us know we are making a second attempt
                    makeBam2(data, bam.bai, indexChunks, callback, true);
                } else {
                    // We've attempted x.bam.bai & x.bai and nothing worked
                    callback(null, result);
                }
            } else {
                bam.data.slice(0, minBlockIndex).fetch(parseBamHeader);
            }
        }); // Timeout on first request to catch Chrome mixed-content error.
    } else {
        var chunks = bam.indexChunks.chunks;
        bam.indices = [];
        for (var i = 0; i < chunks.length; i++) {
            bam.indices[i] = null; // To be filled out lazily as needed
        }
        bam.data.slice(0, minBlockIndex).fetch(parseBamHeader);
    }
}

BamFile.prototype.blocksForRange = function (refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [],
        otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p + 4);
        //        dlog('bin=' + bin + '; nchnk=' + nchnk);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p += nchnk * 16;
        }
    }
    // console.log('leafChunks = ' + JSON.stringify(leafChunks));
    // console.log('otherChunks = ' + JSON.stringify(otherChunks));

    var nintv = readInt(index, p);
    // console.log('nintv=' + nintv);
    var lowest = null;
    var minLin = Math.min(min >> 14, nintv - 1),
        maxLin = Math.min(max >> 14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb = readVob(index, p + 4 + i * 8);
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || lb.block == lowest.block && lb.offset < lowest.offset) {
            lowest = lb;
        }
    }
    // console.log('Lowest LB = ' + lowest);

    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block > lowest.block || chnk.maxv.block == lowest.block && chnk.maxv.offset >= lowest.offset) {
                prunedOtherChunks.push(chnk);
            }
        }
    }
    // console.log('prunedOtherChunks = ' + JSON.stringify(prunedOtherChunks));
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function (c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) {
                    // no point splitting mid-block
                    cur = new Chunk(cur.minv, nc.maxv);
                } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }
    // console.log('mergedChunks = ' + JSON.stringify(mergedChunks));

    return mergedChunks;
};

BamFile.prototype.fetch = function (chr, min, max, callback, opts) {
    var thisB = this;
    opts = opts || {};

    var chrId = this.chrToIndex[chr];
    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        // Fetch this portion of the BAI if it hasn't been loaded yet.
        if (this.indices[chrId] === null && this.indexChunks.chunks[chrId]) {
            var start_stop = this.indexChunks.chunks[chrId];
            return this.bai.slice(start_stop[0], start_stop[1]).fetch(function (data) {
                var buffer = new Uint8Array(data);
                this.indices[chrId] = buffer;
                return this.fetch(chr, min, max, callback, opts);
            }.bind(this));
        }

        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }

    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1 << 16); // *sigh*
            // console.log('fetching ' + fetchMin + ':' + fetchMax);
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function (r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            var finished = thisB.readBamRecords(ba, chunks[index].minv.offset, records, min, max, chrId, opts);
            data = null;
            ++index;
            if (finished) return callback(records);else return tramp();
        }
    }
    tramp();
};

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

function BamRecord() {}

BamFile.prototype.readBamRecords = function (ba, offset, sink, min, max, chrId, opts) {
    while (true) {
        var blockSize = readInt(ba, offset);
        var blockEnd = offset + blockSize + 4;
        if (blockEnd > ba.length) {
            return false;
        }

        var record = new BamRecord();

        var refID = readInt(ba, offset + 4);
        var pos = readInt(ba, offset + 8);

        var bmn = readInt(ba, offset + 12);
        var bin = (bmn & 0xffff0000) >> 16;
        var mq = (bmn & 0xff00) >> 8;
        var nl = bmn & 0xff;

        var flag_nc = readInt(ba, offset + 16);
        var flag = (flag_nc & 0xffff0000) >> 16;
        var nc = flag_nc & 0xffff;

        var lseq = readInt(ba, offset + 20);

        var nextRef = readInt(ba, offset + 24);
        var nextPos = readInt(ba, offset + 28);

        var tlen = readInt(ba, offset + 32);

        record.segment = this.indexToChr[refID];
        record.flag = flag;
        record.pos = pos;
        record.mq = mq;
        if (opts.light) record.seqLength = lseq;

        if (!opts.light || opts.includeName) {
            var readName = '';
            for (var j = 0; j < nl - 1; ++j) {
                readName += String.fromCharCode(ba[offset + 36 + j]);
            }
            record.readName = readName;
        }

        if (!opts.light) {
            if (nextRef >= 0) {
                record.nextSegment = this.indexToChr[nextRef];
                record.nextPos = nextPos;
            }

            var p = offset + 36 + nl;

            var cigar = '';
            for (var c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                cigar = cigar + (cigop >> 4) + CIGAR_DECODER[cigop & 0xf];
                p += 4;
            }
            record.cigar = cigar;

            var seq = '';
            var seqBytes = lseq + 1 >> 1;
            for (var j = 0; j < seqBytes; ++j) {
                var sb = ba[p + j];
                seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
                if (seq.length < lseq) seq += SEQRET_DECODER[sb & 0x0f];
            }
            p += seqBytes;
            record.seq = seq;

            var qseq = '';
            for (var j = 0; j < lseq; ++j) {
                qseq += String.fromCharCode(ba[p + j] + 33);
            }
            p += lseq;
            record.quals = qseq;

            while (p < blockEnd) {
                var tag = String.fromCharCode(ba[p], ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type == 'i' || type == 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type == 'c' || type == 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type == 's' || type == 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type == 'f') {
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type == 'Z' || type == 'H') {
                    p += 3;
                    value = '';
                    for (;;) {
                        var cc = ba[p++];
                        if (cc == 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else if (type == 'B') {
                    var atype = String.fromCharCode(ba[p + 3]);
                    var alen = readInt(ba, p + 4);
                    var elen;
                    var reader;
                    if (atype == 'i' || atype == 'I' || atype == 'f') {
                        elen = 4;
                        if (atype == 'f') reader = readFloat;else reader = readInt;
                    } else if (atype == 's' || atype == 'S') {
                        elen = 2;
                        reader = readShort;
                    } else if (atype == 'c' || atype == 'C') {
                        elen = 1;
                        reader = readByte;
                    } else {
                        throw 'Unknown array type ' + atype;
                    }

                    p += 8;
                    value = [];
                    for (var i = 0; i < alen; ++i) {
                        value.push(reader(ba, p));
                        p += elen;
                    }
                } else {
                    throw 'Unknown type ' + type;
                }
                record[tag] = value;
            }
        }

        if (!min || record.pos <= max && record.pos + lseq >= min) {
            if (chrId === undefined || refID == chrId) {
                sink.push(record);
            }
        }
        if (record.pos > max) {
            return true;
        }
        offset = blockEnd;
    }

    // Exits via top of loop.
};

if (typeof module !== 'undefined') {
    module.exports = {
        makeBam: makeBam,
        BAM_MAGIC: BAM_MAGIC,
        BAI_MAGIC: BAI_MAGIC,
        BamFlags: BamFlags
    };
}

},{"./bin":3,"./lh3utils":8,"./spans":10}],2:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// bigwig.js: indexed binary WIG (and BED) files
//

"use strict";

if (typeof require !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var das = require('./das');
    var DASFeature = das.DASFeature;
    var DASGroup = das.DASGroup;

    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var bin = require('./bin');
    var readInt = bin.readInt;

    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

var BIG_WIG_MAGIC = 0x888FFC26;
var BIG_WIG_MAGIC_BE = 0x26FC8F88;
var BIG_BED_MAGIC = 0x8789F2EB;
var BIG_BED_MAGIC_BE = 0xEBF28987;

var BIG_WIG_TYPE_GRAPH = 1;
var BIG_WIG_TYPE_VSTEP = 2;
var BIG_WIG_TYPE_FSTEP = 3;

var M1 = 256;
var M2 = 256 * 256;
var M3 = 256 * 256 * 256;
var M4 = 256 * 256 * 256 * 256;

var BED_COLOR_REGEXP = new RegExp("^[0-9]+,[0-9]+,[0-9]+");

function bwg_readOffset(ba, o) {
    var offset = ba[o] + ba[o + 1] * M1 + ba[o + 2] * M2 + ba[o + 3] * M3 + ba[o + 4] * M4;
    return offset;
}

function BigWig() {}

BigWig.prototype.readChromTree = function (callback) {
    var thisB = this;
    this.chromsToIDs = {};
    this.idsToChroms = {};
    this.maxID = 0;

    var udo = this.unzoomedDataOffset;
    var eb = udo - this.chromTreeOffset & 3;
    udo = udo + 4 - eb;

    this.data.slice(this.chromTreeOffset, udo - this.chromTreeOffset).fetch(function (bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        var bptReadNode = function bptReadNode(offset) {
            var nodeType = ba[offset];
            var cnt = sa[offset / 2 + 1];
            offset += 4;
            for (var n = 0; n < cnt; ++n) {
                if (nodeType == 0) {
                    offset += keySize;
                    var childOffset = bwg_readOffset(ba, offset);
                    offset += 8;
                    childOffset -= thisB.chromTreeOffset;
                    bptReadNode(childOffset);
                } else {
                    var key = '';
                    for (var ki = 0; ki < keySize; ++ki) {
                        var charCode = ba[offset++];
                        if (charCode != 0) {
                            key += String.fromCharCode(charCode);
                        }
                    }
                    var chromId = ba[offset + 3] << 24 | ba[offset + 2] << 16 | ba[offset + 1] << 8 | ba[offset + 0];
                    var chromSize = ba[offset + 7] << 24 | ba[offset + 6] << 16 | ba[offset + 5] << 8 | ba[offset + 4];
                    offset += 8;

                    thisB.chromsToIDs[key] = chromId;
                    if (key.indexOf('chr') == 0) {
                        thisB.chromsToIDs[key.substr(3)] = chromId;
                    }
                    thisB.idsToChroms[chromId] = key;
                    thisB.maxID = Math.max(thisB.maxID, chromId);
                }
            }
        };
        bptReadNode(rootNodeOffset);

        callback(thisB);
    });
};

function BigWigView(bwg, cirTreeOffset, cirTreeLength, isSummary) {
    this.bwg = bwg;
    this.cirTreeOffset = cirTreeOffset;
    this.cirTreeLength = cirTreeLength;
    this.isSummary = isSummary;
}

BigWigView.prototype.readWigData = function (chrName, min, max, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.readWigDataById(chr, min, max, callback);
    }
};

BigWigView.prototype.readWigDataById = function (chr, min, max, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function (result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.readWigDataById(chr, min, max, callback);
        });
        return;
    }

    var blocksToFetch = [];
    var outstanding = 0;

    var beforeBWG = Date.now();

    var filter = function filter(chromId, fmin, fmax, toks) {
        return (chr < 0 || chromId == chr) && fmin <= max && fmax >= min;
    };

    var cirFobRecur = function cirFobRecur(offset, level) {
        if (thisB.bwg.instrument) console.log('level=' + level + '; offset=' + offset + '; time=' + (Date.now() | 0));

        outstanding += offset.length;

        if (offset.length == 1 && offset[0] - thisB.cirTreeOffset == 48 && thisB.cachedCirRoot) {
            cirFobRecur2(thisB.cachedCirRoot, 0, level);
            --outstanding;
            if (outstanding == 0) {
                thisB.fetchFeatures(filter, blocksToFetch, callback);
            }
            return;
        }

        var maxCirBlockSpan = 4 + thisB.cirBlockSize * 32; // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }

        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    };

    var cirFobStartFetch = function cirFobStartFetch(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function (resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);

                    if (offset[i] - thisB.cirTreeOffset == 48 && offset[i] - fr.min() == 0) thisB.cachedCirRoot = resultBuffer;

                    --outstanding;
                    if (outstanding == 0) {
                        thisB.fetchFeatures(filter, blocksToFetch, callback);
                    }
                }
            }
        });
    };

    var cirFobRecur2 = function cirFobRecur2(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData, 0, cirBlockData.byteLength >> 1);
        var la = new Int32Array(cirBlockData, 0, cirBlockData.byteLength >> 2);

        var isLeaf = ba[offset];
        var cnt = sa[offset / 2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset / 4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset + 16);
                var blockSize = bwg_readOffset(ba, offset + 24);
                if ((chr < 0 || startChrom < chr || startChrom == chr && startBase <= max) && (chr < 0 || endChrom > chr || endChrom == chr && endBase >= min)) {
                    blocksToFetch.push({ offset: blockOffset, size: blockSize });
                }
                offset += 32;
            }
        } else {
            var recurOffsets = [];
            for (var i = 0; i < cnt; ++i) {
                var lo = offset / 4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset + 16);
                if ((chr < 0 || startChrom < chr || startChrom == chr && startBase <= max) && (chr < 0 || endChrom > chr || endChrom == chr && endBase >= min)) {
                    recurOffsets.push(blockOffset);
                }
                offset += 24;
            }
            if (recurOffsets.length > 0) {
                cirFobRecur(recurOffsets, level + 1);
            }
        }
    };

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
};

BigWigView.prototype.fetchFeatures = function (filter, blocksToFetch, callback) {
    var thisB = this;

    blocksToFetch.sort(function (b0, b1) {
        return (b0.offset | 0) - (b1.offset | 0);
    });

    if (blocksToFetch.length == 0) {
        callback([]);
    } else {
        var features = [];
        var createFeature = function createFeature(chr, fmin, fmax, opts) {
            if (!opts) {
                opts = {};
            }

            var f = new DASFeature();
            f._chromId = chr;
            f.segment = thisB.bwg.idsToChroms[chr];
            f.min = fmin;
            f.max = fmax;
            f.type = thisB.bwg.type;

            for (var k in opts) {
                f[k] = opts[k];
            }

            features.push(f);
        };

        var tramp = function tramp() {
            if (blocksToFetch.length == 0) {
                var afterBWG = Date.now();
                // dlog('BWG fetch took ' + (afterBWG - beforeBWG) + 'ms');
                callback(features);
                return; // just in case...
            } else {
                var block = blocksToFetch[0];
                if (block.data) {
                    thisB.parseFeatures(block.data, createFeature, filter);
                    blocksToFetch.splice(0, 1);
                    tramp();
                } else {
                    var fetchStart = block.offset;
                    var fetchSize = block.size;
                    var bi = 1;
                    while (bi < blocksToFetch.length && blocksToFetch[bi].offset == fetchStart + fetchSize) {
                        fetchSize += blocksToFetch[bi].size;
                        ++bi;
                    }

                    thisB.bwg.data.slice(fetchStart, fetchSize).fetch(function (result) {
                        var offset = 0;
                        var bi = 0;
                        while (offset < fetchSize) {
                            var fb = blocksToFetch[bi];

                            var data;
                            if (thisB.bwg.uncompressBufSize > 0) {
                                data = jszlib_inflate_buffer(result, offset + 2, fb.size - 2);
                            } else {
                                var tmp = new Uint8Array(fb.size); // FIXME is this really the best we can do?
                                arrayCopy(new Uint8Array(result, offset, fb.size), 0, tmp, 0, fb.size);
                                data = tmp.buffer;
                            }
                            fb.data = data;

                            offset += fb.size;
                            ++bi;
                        }
                        tramp();
                    });
                }
            }
        };
        tramp();
    }
};

BigWigView.prototype.parseFeatures = function (data, createFeature, filter) {
    var ba = new Uint8Array(data);

    if (this.isSummary) {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var itemCount = data.byteLength / 32;
        for (var i = 0; i < itemCount; ++i) {
            var chromId = la[i * 8];
            var start = la[i * 8 + 1];
            var end = la[i * 8 + 2];
            var validCnt = la[i * 8 + 3];
            var minVal = fa[i * 8 + 4];
            var maxVal = fa[i * 8 + 5];
            var sumData = fa[i * 8 + 6];
            var sumSqData = fa[i * 8 + 7];

            if (filter(chromId, start + 1, end)) {
                var summaryOpts = { type: 'bigwig', score: sumData / validCnt, maxScore: maxVal };
                if (this.bwg.type == 'bigbed') {
                    summaryOpts.type = 'density';
                }
                createFeature(chromId, start + 1, end, summaryOpts);
            }
        }
    } else if (this.bwg.type == 'bigwig') {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var chromId = la[0];
        var blockStart = la[1];
        var blockEnd = la[2];
        var itemStep = la[3];
        var itemSpan = la[4];
        var blockType = ba[20];
        var itemCount = sa[11];

        if (blockType == BIG_WIG_TYPE_FSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var score = fa[i + 6];
                var fmin = blockStart + i * itemStep + 1,
                    fmax = blockStart + i * itemStep + itemSpan;
                if (filter(chromId, fmin, fmax)) createFeature(chromId, fmin, fmax, { score: score });
            }
        } else if (blockType == BIG_WIG_TYPE_VSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[i * 2 + 6] + 1;
                var end = start + itemSpan - 1;
                var score = fa[i * 2 + 7];
                if (filter(chromId, start, end)) createFeature(chromId, start, end, { score: score });
            }
        } else if (blockType == BIG_WIG_TYPE_GRAPH) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[i * 3 + 6] + 1;
                var end = la[i * 3 + 7];
                var score = fa[i * 3 + 8];
                if (start > end) {
                    start = end;
                }
                if (filter(chromId, start, end)) createFeature(chromId, start, end, { score: score });
            }
        } else {
            console.log('Currently not handling bwgType=' + blockType);
        }
    } else if (this.bwg.type == 'bigbed') {
        var offset = 0;
        var dfc = this.bwg.definedFieldCount;
        var schema = this.bwg.schema;

        while (offset < ba.length) {
            var chromId = ba[offset + 3] << 24 | ba[offset + 2] << 16 | ba[offset + 1] << 8 | ba[offset + 0];
            var start = ba[offset + 7] << 24 | ba[offset + 6] << 16 | ba[offset + 5] << 8 | ba[offset + 4];
            var end = ba[offset + 11] << 24 | ba[offset + 10] << 16 | ba[offset + 9] << 8 | ba[offset + 8];
            offset += 12;
            var rest = '';
            while (true) {
                var ch = ba[offset++];
                if (ch != 0) {
                    rest += String.fromCharCode(ch);
                } else {
                    break;
                }
            }

            var featureOpts = {};

            var bedColumns;
            if (rest.length > 0) {
                bedColumns = rest.split('\t');
            } else {
                bedColumns = [];
            }
            if (bedColumns.length > 0 && dfc > 3) {
                featureOpts.label = bedColumns[0];
            }
            if (bedColumns.length > 1 && dfc > 4) {
                var score = parseInt(bedColumns[1]);
                if (!isNaN(score)) featureOpts.score = score;
            }
            if (bedColumns.length > 2 && dfc > 5) {
                featureOpts.orientation = bedColumns[2];
            }
            if (bedColumns.length > 5 && dfc > 8) {
                var color = bedColumns[5];
                if (BED_COLOR_REGEXP.test(color)) {
                    featureOpts.itemRgb = 'rgb(' + color + ')';
                }
            }

            if (bedColumns.length > dfc - 3 && schema) {
                for (var col = dfc - 3; col < bedColumns.length; ++col) {
                    featureOpts[schema.fields[col + 3].name] = bedColumns[col];
                }
            }

            if (filter(chromId, start + 1, end, bedColumns)) {
                if (dfc < 12) {
                    createFeature(chromId, start + 1, end, featureOpts);
                } else {
                    var thickStart = bedColumns[3] | 0;
                    var thickEnd = bedColumns[4] | 0;
                    var blockCount = bedColumns[6] | 0;
                    var blockSizes = bedColumns[7].split(',');
                    var blockStarts = bedColumns[8].split(',');

                    if (featureOpts.exonFrames) {
                        var exonFrames = featureOpts.exonFrames.split(',');
                        featureOpts.exonFrames = undefined;
                    }

                    featureOpts.type = 'transcript';
                    var grp = new DASGroup();
                    for (var k in featureOpts) {
                        grp[k] = featureOpts[k];
                    }
                    grp.id = bedColumns[0];
                    grp.segment = this.bwg.idsToChroms[chromId];
                    grp.min = start + 1;
                    grp.max = end;
                    grp.notes = [];
                    featureOpts.groups = [grp];

                    // Moving towards using bigGenePred model, but will
                    // still support old Dalliance-style BED12+gene-name for the
                    // foreseeable future.
                    if (bedColumns.length > 9) {
                        var geneId = featureOpts.geneName || bedColumns[9];
                        var geneName = geneId;
                        if (bedColumns.length > 10) {
                            geneName = bedColumns[10];
                        }
                        if (featureOpts.geneName2) geneName = featureOpts.geneName2;

                        var gg = shallowCopy(grp);
                        gg.id = geneId;
                        gg.label = geneName;
                        gg.type = 'gene';
                        featureOpts.groups.push(gg);
                    }

                    var spanList = [];
                    for (var b = 0; b < blockCount; ++b) {
                        var bmin = (blockStarts[b] | 0) + start;
                        var bmax = bmin + (blockSizes[b] | 0);
                        var span = new Range(bmin, bmax);
                        spanList.push(span);
                    }
                    var spans = union(spanList);

                    var tsList = spans.ranges();
                    for (var s = 0; s < tsList.length; ++s) {
                        var ts = tsList[s];
                        createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                    }

                    if (thickEnd > thickStart) {
                        var codingRegion = featureOpts.orientation == '+' ? new Range(thickStart, thickEnd + 3) : new Range(thickStart - 3, thickEnd);
                        // +/- 3 to account for stop codon

                        var tl = intersection(spans, codingRegion);
                        if (tl) {
                            featureOpts.type = 'translation';
                            var tlList = tl.ranges();
                            var readingFrame = 0;

                            var tlOffset = 0;
                            while (tlList[0].min() > tsList[tlOffset].max()) {
                                tlOffset++;
                            }for (var s = 0; s < tlList.length; ++s) {
                                // Record reading frame for every exon
                                var index = s;
                                if (featureOpts.orientation == '-') index = tlList.length - s - 1;
                                var ts = tlList[index];
                                featureOpts.readframe = readingFrame;
                                if (exonFrames) {
                                    var brf = parseInt(exonFrames[index + tlOffset]);
                                    if (typeof brf === 'number' && brf >= 0 && brf <= 2) {
                                        featureOpts.readframe = brf;
                                        featureOpts.readframeExplicit = true;
                                    }
                                }
                                var length = ts.max() - ts.min();
                                readingFrame = (readingFrame + length) % 3;
                                createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                            }
                        }
                    }
                }
            }
        }
    } else {
        throw Error("Don't know what to do with " + this.bwg.type);
    }
};

//
// nasty cut/paste, should roll back in!
//

BigWigView.prototype.getFirstAdjacent = function (chrName, pos, dir, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.getFirstAdjacentById(chr, pos, dir, callback);
    }
};

BigWigView.prototype.getFirstAdjacentById = function (chr, pos, dir, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function (result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.getFirstAdjacentById(chr, pos, dir, callback);
        });
        return;
    }

    var blockToFetch = null;
    var bestBlockChr = -1;
    var bestBlockOffset = -1;

    var outstanding = 0;

    var beforeBWG = Date.now();

    var cirFobRecur = function cirFobRecur(offset, level) {
        outstanding += offset.length;

        var maxCirBlockSpan = 4 + thisB.cirBlockSize * 32; // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }

        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    };

    var cirFobStartFetch = function cirFobStartFetch(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function (resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);
                    --outstanding;
                    if (outstanding == 0) {
                        if (!blockToFetch) {
                            if (dir > 0 && (chr != 0 || pos > 0)) {
                                return thisB.getFirstAdjacentById(0, 0, dir, callback);
                            } else if (dir < 0 && (chr != thisB.bwg.maxID || pos < 1000000000)) {
                                return thisB.getFirstAdjacentById(thisB.bwg.maxID, 1000000000, dir, callback);
                            }
                            return callback([]);
                        }

                        thisB.fetchFeatures(function (chrx, fmin, fmax, toks) {
                            return dir < 0 && (chrx < chr || fmax < pos) || dir > 0 && (chrx > chr || fmin > pos);
                        }, [blockToFetch], function (features) {
                            var bestFeature = null;
                            var bestChr = -1;
                            var bestPos = -1;
                            for (var fi = 0; fi < features.length; ++fi) {
                                var f = features[fi];
                                var chrx = f._chromId,
                                    fmin = f.min,
                                    fmax = f.max;
                                if (bestFeature == null || dir < 0 && (chrx > bestChr || fmax > bestPos) || dir > 0 && (chrx < bestChr || fmin < bestPos)) {
                                    bestFeature = f;
                                    bestPos = dir < 0 ? fmax : fmin;
                                    bestChr = chrx;
                                }
                            }

                            if (bestFeature != null) return callback([bestFeature]);else return callback([]);
                        });
                    }
                }
            }
        });
    };

    var cirFobRecur2 = function cirFobRecur2(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset / 2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset / 4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset + 16);
                var blockSize = bwg_readOffset(ba, offset + 24);
                if (dir < 0 && (startChrom < chr || startChrom == chr && startBase <= pos) || dir > 0 && (endChrom > chr || endChrom == chr && endBase >= pos)) {
                    // console.log('Got an interesting block: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                    if (/_random/.exec(thisB.bwg.idsToChroms[startChrom])) {
                        // dlog('skipping random: ' + thisB.bwg.idsToChroms[startChrom]);
                    } else if (blockToFetch == null || dir < 0 && (endChrom > bestBlockChr || endChrom == bestBlockChr && endBase > bestBlockOffset) || dir > 0 && (startChrom < bestBlockChr || startChrom == bestBlockChr && startBase < bestBlockOffset)) {
                        //                        dlog('best is: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                        blockToFetch = { offset: blockOffset, size: blockSize };
                        bestBlockOffset = dir < 0 ? endBase : startBase;
                        bestBlockChr = dir < 0 ? endChrom : startChrom;
                    }
                }
                offset += 32;
            }
        } else {
            var bestRecur = -1;
            var bestPos = -1;
            var bestChr = -1;
            for (var i = 0; i < cnt; ++i) {
                var lo = offset / 4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = la[lo + 4] << 32 | la[lo + 5];
                if (dir < 0 && (startChrom < chr || startChrom == chr && startBase <= pos) && endChrom >= chr || dir > 0 && (endChrom > chr || endChrom == chr && endBase >= pos) && startChrom <= chr) {
                    if (bestRecur < 0 || endBase > bestPos) {
                        bestRecur = blockOffset;
                        bestPos = dir < 0 ? endBase : startBase;
                        bestChr = dir < 0 ? endChrom : startChrom;
                    }
                }
                offset += 24;
            }
            if (bestRecur >= 0) {
                cirFobRecur([bestRecur], level + 1);
            }
        }
    };

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
};

BigWig.prototype.readWigData = function (chrName, min, max, callback) {
    this.getUnzoomedView().readWigData(chrName, min, max, callback);
};

BigWig.prototype.getUnzoomedView = function () {
    if (!this.unzoomedView) {
        var cirLen = 4000;
        var nzl = this.zoomLevels[0];
        if (nzl) {
            cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
        }
        this.unzoomedView = new BigWigView(this, this.unzoomedIndexOffset, cirLen, false);
    }
    return this.unzoomedView;
};

BigWig.prototype.getZoomedView = function (z) {
    var zh = this.zoomLevels[z];
    if (!zh.view) {
        zh.view = new BigWigView(this, zh.indexOffset, /* this.zoomLevels[z + 1].dataOffset - zh.indexOffset */4000, true);
    }
    return zh.view;
};

function makeBwg(data, callback, name) {
    var bwg = new BigWig();
    bwg.data = data;
    bwg.name = name;
    bwg.data.slice(0, 512).salted().fetch(function (result) {
        if (!result) {
            return callback(null, "Couldn't fetch file");
        }

        var header = result;
        var ba = new Uint8Array(header);
        var sa = new Int16Array(header);
        var la = new Int32Array(header);
        var magic = ba[0] + M1 * ba[1] + M2 * ba[2] + M3 * ba[3];
        if (magic == BIG_WIG_MAGIC) {
            bwg.type = 'bigwig';
        } else if (magic == BIG_BED_MAGIC) {
            bwg.type = 'bigbed';
        } else if (magic == BIG_WIG_MAGIC_BE || magic == BIG_BED_MAGIC_BE) {
            return callback(null, "Currently don't support big-endian BBI files");
        } else {
            return callback(null, "Not a supported format, magic=0x" + magic.toString(16));
        }

        bwg.version = sa[2]; // 4
        bwg.numZoomLevels = sa[3]; // 6
        bwg.chromTreeOffset = bwg_readOffset(ba, 8);
        bwg.unzoomedDataOffset = bwg_readOffset(ba, 16);
        bwg.unzoomedIndexOffset = bwg_readOffset(ba, 24);
        bwg.fieldCount = sa[16]; // 32
        bwg.definedFieldCount = sa[17]; // 34
        bwg.asOffset = bwg_readOffset(ba, 36);
        bwg.totalSummaryOffset = bwg_readOffset(ba, 44);
        bwg.uncompressBufSize = la[13]; // 52
        bwg.extHeaderOffset = bwg_readOffset(ba, 56);

        bwg.zoomLevels = [];
        for (var zl = 0; zl < bwg.numZoomLevels; ++zl) {
            var zlReduction = la[zl * 6 + 16];
            var zlData = bwg_readOffset(ba, zl * 24 + 72);
            var zlIndex = bwg_readOffset(ba, zl * 24 + 80);
            bwg.zoomLevels.push({ reduction: zlReduction, dataOffset: zlData, indexOffset: zlIndex });
        }

        bwg.readChromTree(function () {
            bwg.getAutoSQL(function (as) {
                bwg.schema = as;
                return callback(bwg);
            });
        });
    }, { timeout: 5000 }); // Potential timeout on first request to catch mixed-content errors on
    // Chromium.
}

BigWig.prototype._tsFetch = function (zoom, chr, min, max, callback) {
    var bwg = this;
    if (zoom >= this.zoomLevels.length - 1) {
        if (!this.topLevelReductionCache) {
            this.getZoomedView(this.zoomLevels.length - 1).readWigDataById(-1, 0, 300000000, function (feats) {
                bwg.topLevelReductionCache = feats;
                return bwg._tsFetch(zoom, chr, min, max, callback);
            });
        } else {
            var f = [];
            var c = this.topLevelReductionCache;
            for (var fi = 0; fi < c.length; ++fi) {
                if (c[fi]._chromId == chr) {
                    f.push(c[fi]);
                }
            }
            return callback(f);
        }
    } else {
        var view;
        if (zoom < 0) {
            view = this.getUnzoomedView();
        } else {
            view = this.getZoomedView(zoom);
        }
        return view.readWigDataById(chr, min, max, callback);
    }
};

BigWig.prototype.thresholdSearch = function (chrName, referencePoint, dir, threshold, callback) {
    dir = dir < 0 ? -1 : 1;
    var bwg = this;
    var initialChr = this.chromsToIDs[chrName];
    var candidates = [{ chrOrd: 0, chr: initialChr, zoom: bwg.zoomLevels.length - 4, min: 0, max: 300000000, fromRef: true }];
    for (var i = 1; i <= this.maxID + 1; ++i) {
        var chrId = (initialChr + dir * i) % (this.maxID + 1);
        if (chrId < 0) chrId += this.maxID + 1;
        candidates.push({ chrOrd: i, chr: chrId, zoom: bwg.zoomLevels.length - 1, min: 0, max: 300000000 });
    }

    function fbThresholdSearchRecur() {
        if (candidates.length == 0) {
            return callback(null);
        }
        candidates.sort(function (c1, c2) {
            var d = c1.zoom - c2.zoom;
            if (d != 0) return d;

            d = c1.chrOrd - c2.chrOrd;
            if (d != 0) return d;else return c1.min - c2.min * dir;
        });

        var candidate = candidates.splice(0, 1)[0];
        bwg._tsFetch(candidate.zoom, candidate.chr, candidate.min, candidate.max, function (feats) {
            var rp = dir > 0 ? 0 : 300000000;
            if (candidate.fromRef) rp = referencePoint;

            for (var fi = 0; fi < feats.length; ++fi) {
                var f = feats[fi];
                var score;
                if (f.maxScore != undefined) score = f.maxScore;else score = f.score;

                if (dir > 0) {
                    if (score > threshold) {
                        if (candidate.zoom < 0) {
                            if (f.min > rp) return callback(f);
                        } else if (f.max > rp) {
                            candidates.push({ chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef });
                        }
                    }
                } else {
                    if (score > threshold) {
                        if (candidate.zoom < 0) {
                            if (f.max < rp) return callback(f);
                        } else if (f.min < rp) {
                            candidates.push({ chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef });
                        }
                    }
                }
            }
            fbThresholdSearchRecur();
        });
    }

    fbThresholdSearchRecur();
};

BigWig.prototype.getAutoSQL = function (callback) {
    var thisB = this;
    if (!this.asOffset) return callback(null);

    this.data.slice(this.asOffset, 2048).fetch(function (result) {
        var ba = new Uint8Array(result);
        var s = '';
        for (var i = 0; i < ba.length; ++i) {
            if (ba[i] == 0) break;
            s += String.fromCharCode(ba[i]);
        }

        /* 
         * Quick'n'dirty attempt to parse autoSql format.
         * See: http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/059/5949/5949l2.html
         */

        var header_re = /(\w+)\s+(\w+)\s+("([^"]+)")?\s+\(\s*/;
        var field_re = /([\w\[\]]+)\s+(\w+)\s*;\s*("([^"]+)")?\s*/g;

        var headerMatch = header_re.exec(s);
        if (headerMatch) {
            var as = {
                declType: headerMatch[1],
                name: headerMatch[2],
                comment: headerMatch[4],

                fields: []
            };

            s = s.substring(headerMatch[0]);
            for (var m = field_re.exec(s); m != null; m = field_re.exec(s)) {
                as.fields.push({ type: m[1],
                    name: m[2],
                    comment: m[4] });
            }

            return callback(as);
        }
    });
};

BigWig.prototype.getExtraIndices = function (callback) {
    var thisB = this;
    if (this.version < 4 || this.extHeaderOffset == 0 || this.type != 'bigbed') {
        return callback(null);
    } else {
        this.data.slice(this.extHeaderOffset, 64).fetch(function (result) {
            if (!result) {
                return callback(null, "Couldn't fetch extension header");
            }

            var ba = new Uint8Array(result);
            var sa = new Int16Array(result);
            var la = new Int32Array(result);

            var extHeaderSize = sa[0];
            var extraIndexCount = sa[1];
            var extraIndexListOffset = bwg_readOffset(ba, 4);

            if (extraIndexCount == 0) {
                return callback(null);
            }

            // FIXME 20byte records only make sense for single-field indices.
            // Right now, these seem to be the only things around, but the format
            // is actually more general.
            thisB.data.slice(extraIndexListOffset, extraIndexCount * 20).fetch(function (eil) {
                if (!eil) {
                    return callback(null, "Couldn't fetch index info");
                }

                var ba = new Uint8Array(eil);
                var sa = new Int16Array(eil);
                var la = new Int32Array(eil);

                var indices = [];
                for (var ii = 0; ii < extraIndexCount; ++ii) {
                    var eiType = sa[ii * 10];
                    var eiFieldCount = sa[ii * 10 + 1];
                    var eiOffset = bwg_readOffset(ba, ii * 20 + 4);
                    var eiField = sa[ii * 10 + 8];
                    var index = new BBIExtraIndex(thisB, eiType, eiFieldCount, eiOffset, eiField);
                    indices.push(index);
                }
                callback(indices);
            });
        });
    }
};

function BBIExtraIndex(bbi, type, fieldCount, offset, field) {
    this.bbi = bbi;
    this.type = type;
    this.fieldCount = fieldCount;
    this.offset = offset;
    this.field = field;
}

BBIExtraIndex.prototype.lookup = function (name, callback) {
    var thisB = this;

    this.bbi.data.slice(this.offset, 32).fetch(function (bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        function bptReadNode(nodeOffset) {
            thisB.bbi.data.slice(nodeOffset, 4 + blockSize * (keySize + valSize)).fetch(function (node) {
                var ba = new Uint8Array(node);
                var sa = new Uint16Array(node);
                var la = new Uint32Array(node);

                var nodeType = ba[0];
                var cnt = sa[1];

                var offset = 4;
                if (nodeType == 0) {
                    var lastChildOffset = null;
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }

                        var childOffset = bwg_readOffset(ba, offset);
                        offset += 8;

                        if (name.localeCompare(key) < 0 && lastChildOffset) {
                            bptReadNode(lastChildOffset);
                            return;
                        }
                        lastChildOffset = childOffset;
                    }
                    bptReadNode(lastChildOffset);
                } else {
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }

                        // Specific for EI case.
                        if (key == name) {
                            var start = bwg_readOffset(ba, offset);
                            var length = readInt(ba, offset + 8);

                            return thisB.bbi.getUnzoomedView().fetchFeatures(function (chr, min, max, toks) {
                                if (toks && toks.length > thisB.field - 3) return toks[thisB.field - 3] == name;
                            }, [{ offset: start, size: length }], callback);
                        }
                        offset += valSize;
                    }
                    return callback([]);
                }
            });
        }

        bptReadNode(thisB.offset + rootNodeOffset);
    });
};

if (typeof module !== 'undefined') {
    module.exports = {
        makeBwg: makeBwg,
        BIG_BED_MAGIC: BIG_BED_MAGIC,
        BIG_WIG_MAGIC: BIG_WIG_MAGIC
    };
}

},{"./bin":3,"./das":5,"./spans":10,"./utils":11,"jszlib":13}],3:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bin.js general binary data support
//

"use strict";

var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

if (typeof require !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;

    var Promise = require('es6-promise').Promise;
}

function BlobFetchable(b) {
    this.blob = b;
}

BlobFetchable.prototype.slice = function (start, length) {
    var b;

    if (this.blob.slice) {
        if (length) {
            b = this.blob.slice(start, start + length);
        } else {
            b = this.blob.slice(start);
        }
    } else {
        if (length) {
            b = this.blob.webkitSlice(start, start + length);
        } else {
            b = this.blob.webkitSlice(start);
        }
    }
    return new BlobFetchable(b);
};

BlobFetchable.prototype.salted = function () {
    return this;
};

if (typeof FileReader !== 'undefined') {
    // console.log('defining async BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function (callback) {
        var reader = new FileReader();
        reader.onloadend = function (ev) {
            callback(bstringToBuffer(reader.result));
        };
        reader.readAsBinaryString(this.blob);
    };
} else {
    // if (console && console.log)
    //    console.log('defining sync BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function (callback) {
        var reader = new FileReaderSync();
        try {
            var res = reader.readAsArrayBuffer(this.blob);
            callback(res);
        } catch (e) {
            callback(null, e);
        }
    };
}

function URLFetchable(url, start, end, opts) {
    if (!opts) {
        if ((typeof start === 'undefined' ? 'undefined' : _typeof(start)) === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = url;
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}

URLFetchable.prototype.slice = function (s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start,
        ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new URLFetchable(this.url, ns, ne, this.opts);
};

var seed = 0;
var isSafari = typeof navigator !== 'undefined' && navigator.userAgent.indexOf('Safari') >= 0 && navigator.userAgent.indexOf('Chrome') < 0;

URLFetchable.prototype.fetchAsText = function (callback) {
    var thisB = this;

    this.getURL().then(function (url) {
        try {
            var req = new XMLHttpRequest();
            var length;
            if ((isSafari || thisB.opts.salt) && url.indexOf('?') < 0) {
                url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + ++seed);
            }
            req.open('GET', url, true);

            if (thisB.end) {
                if (thisB.end - thisB.start > 100000000) {
                    throw 'Monster fetch!';
                }
                req.setRequestHeader('Range', 'bytes=' + thisB.start + '-' + thisB.end);
                length = thisB.end - thisB.start + 1;
            }

            req.onreadystatechange = function () {
                if (req.readyState == 4) {
                    if (req.status == 200 || req.status == 206) {
                        return callback(req.responseText);
                    } else {
                        return callback(null);
                    }
                }
            };
            if (thisB.opts.credentials) {
                req.withCredentials = true;
            }
            req.send();
        } catch (e) {
            return callback(null);
        }
    }).catch(function (err) {
        console.log(err);
        return callback(null, err);
    });
};

URLFetchable.prototype.salted = function () {
    var o = shallowCopy(this.opts);
    o.salt = true;
    return new URLFetchable(this.url, this.start, this.end, o);
};

URLFetchable.prototype.getURL = function () {
    if (this.opts.resolver) {
        return this.opts.resolver(this.url).then(function (urlOrObj) {
            if (typeof urlOrObj === 'string') {
                return urlOrObj;
            } else {
                return urlOrObj.url;
            }
        });
    } else {
        return Promise.resolve(this.url);
    }
};

URLFetchable.prototype.fetch = function (callback, opts) {
    var thisB = this;

    opts = opts || {};
    var attempt = opts.attempt || 1;
    var truncatedLength = opts.truncatedLength;
    if (attempt > 3) {
        return callback(null);
    }

    this.getURL().then(function (url) {
        try {
            var timeout;
            if (opts.timeout && !thisB.opts.credentials) {
                timeout = setTimeout(function () {
                    console.log('timing out ' + url);
                    req.abort();
                    return callback(null, 'Timeout');
                }, opts.timeout);
            }

            var req = new XMLHttpRequest();
            var length;
            if ((isSafari || thisB.opts.salt) && url.indexOf('?') < 0) {
                url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + ++seed);
            }
            req.open('GET', url, true);
            req.overrideMimeType('text/plain; charset=x-user-defined');
            if (thisB.end) {
                if (thisB.end - thisB.start > 100000000) {
                    throw 'Monster fetch!';
                }
                req.setRequestHeader('Range', 'bytes=' + thisB.start + '-' + thisB.end);
                length = thisB.end - thisB.start + 1;
            }
            req.responseType = 'arraybuffer';
            req.onreadystatechange = function () {
                if (req.readyState == 4) {
                    if (timeout) clearTimeout(timeout);
                    if (req.status == 200 || req.status == 206) {
                        if (req.response) {
                            var bl = req.response.byteLength;
                            if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                                return thisB.fetch(callback, { attempt: attempt + 1, truncatedLength: bl });
                            } else {
                                return callback(req.response);
                            }
                        } else if (req.mozResponseArrayBuffer) {
                            return callback(req.mozResponseArrayBuffer);
                        } else {
                            var r = req.responseText;
                            if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                                return thisB.fetch(callback, { attempt: attempt + 1, truncatedLength: r.length });
                            } else {
                                return callback(bstringToBuffer(req.responseText));
                            }
                        }
                    } else {
                        return thisB.fetch(callback, { attempt: attempt + 1 });
                    }
                }
            };
            if (thisB.opts.credentials) {
                req.withCredentials = true;
            }
            req.send();
        } catch (e) {
            return callback(null);
        }
    }).catch(function (err) {
        console.log(err);
        return callback(null, err);
    });
};

function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

// Read from Uint8Array

var convertBuffer = new ArrayBuffer(8);
var ba = new Uint8Array(convertBuffer);
var fa = new Float32Array(convertBuffer);

function readFloat(buf, offset) {
    ba[0] = buf[offset];
    ba[1] = buf[offset + 1];
    ba[2] = buf[offset + 2];
    ba[3] = buf[offset + 3];
    return fa[0];
}

function readInt64(ba, offset) {
    return ba[offset + 7] << 24 | ba[offset + 6] << 16 | ba[offset + 5] << 8 | ba[offset + 4];
}

var M1 = 256,
    M2 = M1 * 256,
    M3 = M2 * 256,
    M4 = M3 * 256,
    M5 = M4 * 256;

function readInt64LE(ba, offset) {
    return ba[offset] + ba[offset + 1] * M1 + ba[offset + 2] * M2 + ba[offset + 3] * M3 + ba[offset + 4] * M4 + ba[offset + 5] * M5;
}

function readInt64BE(ba, offset) {
    return ba[offset + 7] + ba[offset + 6] * M1 + ba[offset + 5] * M2 + ba[offset + 4] * M3 + ba[offset + 3] * M4 + ba[offset + 2] * M5;
}

function readInt(ba, offset) {
    return ba[offset + 3] << 24 | ba[offset + 2] << 16 | ba[offset + 1] << 8 | ba[offset];
}

function readShort(ba, offset) {
    return ba[offset + 1] << 8 | ba[offset];
}

function readByte(ba, offset) {
    return ba[offset];
}

function readIntBE(ba, offset) {
    return ba[offset] << 24 | ba[offset + 1] << 16 | ba[offset + 2] << 8 | ba[offset + 3];
}

// Exports if we are being used as a module

if (typeof module !== 'undefined') {
    module.exports = {
        BlobFetchable: BlobFetchable,
        URLFetchable: URLFetchable,

        readInt: readInt,
        readIntBE: readIntBE,
        readInt64: readInt64,
        readInt64LE: readInt64LE,
        readInt64BE: readInt64BE,
        readShort: readShort,
        readByte: readByte,
        readFloat: readFloat
    };
}

},{"./sha1":9,"./utils":11,"es6-promise":12}],4:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// color.js
//

"use strict";

function DColour(red, green, blue, name) {
    this.red = red | 0;
    this.green = green | 0;
    this.blue = blue | 0;
    if (name) {
        this.name = name;
    }
}

DColour.prototype.toSvgString = function () {
    if (!this.name) {
        this.name = "rgb(" + this.red + "," + this.green + "," + this.blue + ")";
    }

    return this.name;
};

function hex2(x) {
    var y = '00' + x.toString(16);
    return y.substring(y.length - 2);
}

DColour.prototype.toHexString = function () {
    return '#' + hex2(this.red) + hex2(this.green) + hex2(this.blue);
};

var palette = {
    red: new DColour(255, 0, 0, 'red'),
    green: new DColour(0, 255, 0, 'green'),
    blue: new DColour(0, 0, 255, 'blue'),
    yellow: new DColour(255, 255, 0, 'yellow'),
    white: new DColour(255, 255, 255, 'white'),
    black: new DColour(0, 0, 0, 'black'),
    gray: new DColour(180, 180, 180, 'gray'),
    grey: new DColour(180, 180, 180, 'grey'),
    lightskyblue: new DColour(135, 206, 250, 'lightskyblue'),
    lightsalmon: new DColour(255, 160, 122, 'lightsalmon'),
    hotpink: new DColour(255, 105, 180, 'hotpink')
};

var COLOR_RE = new RegExp('^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$');
var CSS_COLOR_RE = /rgb\(([0-9]+),([0-9]+),([0-9]+)\)/;

function dasColourForName(name) {
    var c = palette[name];
    if (!c) {
        var match = COLOR_RE.exec(name);
        if (match) {
            c = new DColour('0x' + match[1] | 0, '0x' + match[2] | 0, '0x' + match[3] | 0, name);
            palette[name] = c;
        } else {
            match = CSS_COLOR_RE.exec(name);
            if (match) {
                c = new DColour(match[1] | 0, match[2] | 0, match[3] | 0, name);
                palette[name] = c;
            } else {
                console.log("couldn't handle color: " + name);
                c = palette.black;
                palette[name] = c;
            }
        }
    }
    return c;
}

function makeColourSteps(steps, stops, colours) {
    var dcolours = [];
    for (var ci = 0; ci < colours.length; ++ci) {
        dcolours.push(dasColourForName(colours[ci]));
    }

    var grad = [];
    STEP_LOOP: for (var si = 0; si < steps; ++si) {
        var rs = 1.0 * si / (steps - 1);
        var score = stops[0] + (stops[stops.length - 1] - stops[0]) * rs;
        for (var i = 0; i < stops.length - 1; ++i) {
            if (score >= stops[i] && score <= stops[i + 1]) {
                var frac = (score - stops[i]) / (stops[i + 1] - stops[i]);
                var ca = dcolours[i];
                var cb = dcolours[i + 1];

                var fill = new DColour(ca.red * (1.0 - frac) + cb.red * frac | 0, ca.green * (1.0 - frac) + cb.green * frac | 0, ca.blue * (1.0 - frac) + cb.blue * frac | 0).toSvgString();
                grad.push(fill);

                continue STEP_LOOP;
            }
        }
        throw 'Bad step';
    }

    return grad;
}

function makeGradient(steps, color1, color2, color3) {
    if (color3) {
        return makeColourSteps(steps, [0, 0.5, 1], [color1, color2, color3]);
    } else {
        return makeColourSteps(steps, [0, 1], [color1, color2]);
    }
}

if (typeof module !== 'undefined') {
    module.exports = {
        makeColourSteps: makeColourSteps,
        makeGradient: makeGradient,
        dasColourForName: dasColourForName
    };
}

},{}],5:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// das.js: queries and low-level data model.
//

"use strict";

if (typeof require !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;
    var pusho = utils.pusho;

    var color = require('./color');
    var makeColourSteps = color.makeColourSteps;
}

var dasLibErrorHandler = function dasLibErrorHandler(errMsg) {
    alert(errMsg);
};
var dasLibRequestQueue = new Array();

function DASSegment(name, start, end, description) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.description = description;
}
DASSegment.prototype.toString = function () {
    return this.name + ':' + this.start + '..' + this.end;
};
DASSegment.prototype.isBounded = function () {
    return this.start && this.end;
};
DASSegment.prototype.toDASQuery = function () {
    var q = 'segment=' + this.name;
    if (this.start && this.end) {
        q += ':' + this.start + ',' + this.end;
    }
    return q;
};

function DASSource(a1, a2) {
    var options;
    if (typeof a1 == 'string') {
        this.uri = a1;
        options = a2 || {};
    } else {
        options = a1 || {};
    }
    for (var k in options) {
        this[k] = options[k];
    }

    if (!this.coords) {
        this.coords = [];
    }
    if (!this.props) {
        this.props = {};
    }

    this.dasBaseURI = this.uri;
    if (this.dasBaseURI && this.dasBaseURI.substr(this.uri.length - 1) != '/') {
        this.dasBaseURI = this.dasBaseURI + '/';
    }
}

DASSource.prototype.getURI = function (uri) {
    if (this.resolver) {
        return this.resolver(uri).then(function (urlOrObj) {
            if (typeof urlOrObj === 'string') {
                return urlOrObj;
            } else {
                return urlOrObj.url;
            }
        });
    } else {
        return Promise.resolve(uri);
    }
};

function DASCoords() {}

function coordsMatch(c1, c2) {
    return c1.taxon == c2.taxon && c1.auth == c2.auth && c1.version == c2.version;
}

//
// DAS 1.6 entry_points command
//

DASSource.prototype.entryPoints = function (callback) {
    var dasURI = this.dasBaseURI + 'entry_points';
    this.doCrossDomainRequest(dasURI, function (responseXML) {
        if (!responseXML) {
            return callback([]);
        }

        var entryPoints = new Array();

        var segs = responseXML.getElementsByTagName('SEGMENT');
        for (var i = 0; i < segs.length; ++i) {
            var seg = segs[i];
            var segId = seg.getAttribute('id');

            var segSize = seg.getAttribute('size');
            var segMin, segMax;
            if (segSize) {
                segMin = 1;segMax = segSize | 0;
            } else {
                segMin = seg.getAttribute('start');
                if (segMin) {
                    segMin |= 0;
                }
                segMax = seg.getAttribute('stop');
                if (segMax) {
                    segMax |= 0;
                }
            }
            var segDesc = null;
            if (seg.firstChild) {
                segDesc = seg.firstChild.nodeValue;
            }
            entryPoints.push(new DASSegment(segId, segMin, segMax, segDesc));
        }
        callback(entryPoints);
    });
};

//
// DAS 1.6 sequence command
// Do we need an option to fall back to the dna command?
//

function DASSequence(name, start, end, alpha, seq) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.alphabet = alpha;
    this.seq = seq;
}

DASSource.prototype.sequence = function (segment, callback) {
    var dasURI = this.dasBaseURI + 'sequence?' + segment.toDASQuery();
    this.doCrossDomainRequest(dasURI, function (responseXML) {
        if (!responseXML) {
            callback([]);
            return;
        } else {
            var seqs = new Array();

            var segs = responseXML.getElementsByTagName('SEQUENCE');
            for (var i = 0; i < segs.length; ++i) {
                var seg = segs[i];
                var segId = seg.getAttribute('id');
                var segMin = seg.getAttribute('start');
                var segMax = seg.getAttribute('stop');
                var segAlpha = 'DNA';
                var segSeq = null;
                if (seg.firstChild) {
                    var rawSeq = seg.firstChild.nodeValue;
                    segSeq = '';
                    var idx = 0;
                    while (true) {
                        var space = rawSeq.indexOf('\n', idx);
                        if (space >= 0) {
                            segSeq += rawSeq.substring(idx, space).toUpperCase();
                            idx = space + 1;
                        } else {
                            segSeq += rawSeq.substring(idx).toUpperCase();
                            break;
                        }
                    }
                }
                seqs.push(new DASSequence(segId, segMin, segMax, segAlpha, segSeq));
            }

            callback(seqs);
        }
    });
};

//
// DAS 1.6 features command
//

function DASFeature() {}

function DASGroup(id) {
    if (id) this.id = id;
}

function DASLink(desc, uri) {
    this.desc = desc;
    this.uri = uri;
}

DASSource.prototype.features = function (segment, options, callback) {
    options = options || {};
    var thisB = this;

    var dasURI;
    if (this.features_uri) {
        dasURI = this.features_uri;
    } else {
        var filters = [];

        if (segment) {
            filters.push(segment.toDASQuery());
        } else if (options.group) {
            var g = options.group;
            if (typeof g == 'string') {
                filters.push('group_id=' + g);
            } else {
                for (var gi = 0; gi < g.length; ++gi) {
                    filters.push('group_id=' + g[gi]);
                }
            }
        }

        if (options.adjacent) {
            var adj = options.adjacent;
            if (typeof adj == 'string') {
                adj = [adj];
            }
            for (var ai = 0; ai < adj.length; ++ai) {
                filters.push('adjacent=' + adj[ai]);
            }
        }

        if (options.type) {
            if (typeof options.type == 'string') {
                filters.push('type=' + options.type);
            } else {
                for (var ti = 0; ti < options.type.length; ++ti) {
                    filters.push('type=' + options.type[ti]);
                }
            }
        }

        if (options.maxbins) {
            filters.push('maxbins=' + options.maxbins);
        }

        if (filters.length > 0) {
            dasURI = this.dasBaseURI + 'features?' + filters.join(';');
        } else {
            callback([], 'No filters specified');
        }
    }

    this.doCrossDomainRequest(dasURI, function (responseXML, req) {
        if (!responseXML) {
            var msg;
            if (req.status == 0) {
                msg = 'server may not support CORS';
            } else {
                msg = 'status=' + req.status;
            }
            callback([], 'Failed request: ' + msg);
            return;
        }
        /*      if (req) {
                    var caps = req.getResponseHeader('X-DAS-Capabilties');
                    if (caps) {
                        alert(caps);
                    }
                } */

        var features = new Array();
        var segmentMap = {};

        var segs = responseXML.getElementsByTagName('SEGMENT');
        for (var si = 0; si < segs.length; ++si) {
            var segmentXML = segs[si];
            var segmentID = segmentXML.getAttribute('id');
            segmentMap[segmentID] = {
                min: segmentXML.getAttribute('start'),
                max: segmentXML.getAttribute('stop')
            };

            var featureXMLs = segmentXML.getElementsByTagName('FEATURE');
            for (var i = 0; i < featureXMLs.length; ++i) {
                var feature = featureXMLs[i];
                var dasFeature = new DASFeature();

                dasFeature.segment = segmentID;
                dasFeature.id = feature.getAttribute('id');
                dasFeature.label = feature.getAttribute('label');

                /*
                                var childNodes = feature.childNodes;
                                for (var c = 0; c < childNodes.length; ++c) {
                                    var cn = childNodes[c];
                                    if (cn.nodeType == Node.ELEMENT_NODE) {
                                        var key = cn.tagName;
                                        //var val = null;
                                        //if (cn.firstChild) {
                                        //   val = cn.firstChild.nodeValue;
                                        //}
                                        dasFeature[key] = 'x';
                                    }
                                } */

                var spos = elementValue(feature, "START");
                var epos = elementValue(feature, "END");
                if ((spos | 0) > (epos | 0)) {
                    dasFeature.min = epos | 0;
                    dasFeature.max = spos | 0;
                } else {
                    dasFeature.min = spos | 0;
                    dasFeature.max = epos | 0;
                }
                {
                    var tec = feature.getElementsByTagName('TYPE');
                    if (tec.length > 0) {
                        var te = tec[0];
                        if (te.firstChild) {
                            dasFeature.type = te.firstChild.nodeValue;
                        }
                        dasFeature.typeId = te.getAttribute('id');
                        dasFeature.typeCv = te.getAttribute('cvId');
                    }
                }
                dasFeature.type = elementValue(feature, "TYPE");
                if (!dasFeature.type && dasFeature.typeId) {
                    dasFeature.type = dasFeature.typeId; // FIXME?
                }

                dasFeature.method = elementValue(feature, "METHOD");
                {
                    var ori = elementValue(feature, "ORIENTATION");
                    if (!ori) {
                        ori = '0';
                    }
                    dasFeature.orientation = ori;
                }
                dasFeature.score = elementValue(feature, "SCORE");
                dasFeature.links = dasLinksOf(feature);
                dasFeature.notes = dasNotesOf(feature);

                var groups = feature.getElementsByTagName("GROUP");
                for (var gi = 0; gi < groups.length; ++gi) {
                    var groupXML = groups[gi];
                    var dasGroup = new DASGroup();
                    dasGroup.type = groupXML.getAttribute('type');
                    dasGroup.id = groupXML.getAttribute('id');
                    dasGroup.links = dasLinksOf(groupXML);
                    dasGroup.notes = dasNotesOf(groupXML);
                    if (!dasFeature.groups) {
                        dasFeature.groups = new Array(dasGroup);
                    } else {
                        dasFeature.groups.push(dasGroup);
                    }
                }

                // Magic notes.  Check with TAD before changing this.
                if (dasFeature.notes) {
                    for (var ni = 0; ni < dasFeature.notes.length; ++ni) {
                        var n = dasFeature.notes[ni];
                        if (n.indexOf('Genename=') == 0) {
                            var gg = new DASGroup();
                            gg.type = 'gene';
                            gg.id = n.substring(9);
                            if (!dasFeature.groups) {
                                dasFeature.groups = new Array(gg);
                            } else {
                                dasFeature.groups.push(gg);
                            }
                        }
                    }
                }

                {
                    var pec = feature.getElementsByTagName('PART');
                    if (pec.length > 0) {
                        var parts = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parts.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parts = parts;
                    }
                }
                {
                    var pec = feature.getElementsByTagName('PARENT');
                    if (pec.length > 0) {
                        var parents = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parents.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parents = parents;
                    }
                }

                features.push(dasFeature);
            }
        }

        callback(features, undefined, segmentMap);
    }, function (err) {
        callback([], err);
    });
};

function DASAlignment(type) {
    this.type = type;
    this.objects = {};
    this.blocks = [];
}

DASSource.prototype.alignments = function (segment, options, callback) {
    var dasURI = this.dasBaseURI + 'alignment?query=' + segment;
    this.doCrossDomainRequest(dasURI, function (responseXML) {
        if (!responseXML) {
            callback([], 'Failed request ' + dasURI);
            return;
        }

        var alignments = [];
        var aliXMLs = responseXML.getElementsByTagName('alignment');
        for (var ai = 0; ai < aliXMLs.length; ++ai) {
            var aliXML = aliXMLs[ai];
            var ali = new DASAlignment(aliXML.getAttribute('alignType'));
            var objXMLs = aliXML.getElementsByTagName('alignObject');
            for (var oi = 0; oi < objXMLs.length; ++oi) {
                var objXML = objXMLs[oi];
                var obj = {
                    id: objXML.getAttribute('intObjectId'),
                    accession: objXML.getAttribute('dbAccessionId'),
                    version: objXML.getAttribute('objectVersion'),
                    dbSource: objXML.getAttribute('dbSource'),
                    dbVersion: objXML.getAttribute('dbVersion')
                };
                ali.objects[obj.id] = obj;
            }

            var blockXMLs = aliXML.getElementsByTagName('block');
            for (var bi = 0; bi < blockXMLs.length; ++bi) {
                var blockXML = blockXMLs[bi];
                var block = {
                    order: blockXML.getAttribute('blockOrder'),
                    segments: []
                };
                var segXMLs = blockXML.getElementsByTagName('segment');
                for (var si = 0; si < segXMLs.length; ++si) {
                    var segXML = segXMLs[si];
                    var seg = {
                        object: segXML.getAttribute('intObjectId'),
                        min: segXML.getAttribute('start'),
                        max: segXML.getAttribute('end'),
                        strand: segXML.getAttribute('strand'),
                        cigar: elementValue(segXML, 'cigar')
                    };
                    block.segments.push(seg);
                }
                ali.blocks.push(block);
            }

            alignments.push(ali);
        }
        callback(alignments);
    });
};

function DASStylesheet() {
    this.styles = [];
}

DASStylesheet.prototype.pushStyle = function (filters, zoom, style) {
    if (!filters) {
        filters = { type: 'default' };
    }
    var styleHolder = shallowCopy(filters);
    if (zoom) {
        styleHolder.zoom = zoom;
    }
    styleHolder.style = style;
    this.styles.push(styleHolder);
};

function DASStyle() {}

function parseGradient(grad) {
    var steps = grad.getAttribute('steps');
    if (steps) {
        steps = steps | 0;
    } else {
        steps = 50;
    }

    var stops = [];
    var colors = [];
    var se = grad.getElementsByTagName('STOP');
    for (var si = 0; si < se.length; ++si) {
        var stop = se[si];
        stops.push(1.0 * stop.getAttribute('score'));
        colors.push(stop.firstChild.nodeValue);
    }

    return makeColourSteps(steps, stops, colors);
}

DASSource.prototype.stylesheet = function (successCB, failureCB) {
    var dasURI,
        creds = this.credentials;
    if (this.stylesheet_uri) {
        dasURI = this.stylesheet_uri;
        creds = false;
    } else {
        dasURI = this.dasBaseURI + 'stylesheet';
    }

    this.getURI(dasURI).then(function (dasURI) {
        doCrossDomainRequest(dasURI, function (responseXML) {
            if (!responseXML) {
                if (failureCB) {
                    failureCB();
                }
                return;
            }
            var stylesheet = new DASStylesheet();
            var typeXMLs = responseXML.getElementsByTagName('TYPE');
            for (var i = 0; i < typeXMLs.length; ++i) {
                var typeStyle = typeXMLs[i];

                var filter = {};
                filter.type = typeStyle.getAttribute('id'); // Am I right in thinking that this makes DASSTYLE XML invalid?  Ugh.
                filter.label = typeStyle.getAttribute('label');
                filter.method = typeStyle.getAttribute('method');
                var glyphXMLs = typeStyle.getElementsByTagName('GLYPH');
                for (var gi = 0; gi < glyphXMLs.length; ++gi) {
                    var glyphXML = glyphXMLs[gi];
                    var zoom = glyphXML.getAttribute('zoom');
                    var glyph = childElementOf(glyphXML);
                    var style = new DASStyle();
                    style.glyph = glyph.localName;
                    var child = glyph.firstChild;

                    while (child) {
                        if (child.nodeType == Node.ELEMENT_NODE) {
                            if (child.localName == 'BGGRAD') {
                                style[child.localName] = parseGradient(child);
                            } else {
                                style[child.localName] = child.firstChild.nodeValue;
                            }
                        }
                        child = child.nextSibling;
                    }
                    stylesheet.pushStyle(filter, zoom, style);
                }
            }
            successCB(stylesheet);
        }, creds);
    }).catch(function (err) {
        console.log(err);
        failureCB();
    });
};

//
// sources command
// 

function DASRegistry(uri, opts) {
    opts = opts || {};
    this.uri = uri;
    this.opts = opts;
}

DASRegistry.prototype.sources = function (callback, failure, opts) {
    if (!opts) {
        opts = {};
    }

    var filters = [];
    if (opts.taxon) {
        filters.push('organism=' + opts.taxon);
    }
    if (opts.auth) {
        filters.push('authority=' + opts.auth);
    }
    if (opts.version) {
        filters.push('version=' + opts.version);
    }
    var quri = this.uri;
    if (filters.length > 0) {
        quri = quri + '?' + filters.join('&'); // '&' as a separator to hack around dasregistry.org bug.
    }

    doCrossDomainRequest(quri, function (responseXML) {
        if (!responseXML && failure) {
            failure();
            return;
        }

        var sources = [];
        var sourceXMLs = responseXML.getElementsByTagName('SOURCE');
        for (var si = 0; si < sourceXMLs.length; ++si) {
            var sourceXML = sourceXMLs[si];
            var versionXMLs = sourceXML.getElementsByTagName('VERSION');
            if (versionXMLs.length < 1) {
                continue;
            }
            var versionXML = versionXMLs[0];

            var coordXMLs = versionXML.getElementsByTagName('COORDINATES');
            var coords = [];
            for (var ci = 0; ci < coordXMLs.length; ++ci) {
                var coordXML = coordXMLs[ci];
                var coord = new DASCoords();
                coord.auth = coordXML.getAttribute('authority');
                coord.taxon = coordXML.getAttribute('taxid');
                coord.version = coordXML.getAttribute('version');
                coords.push(coord);
            }

            var caps = [];
            var capXMLs = versionXML.getElementsByTagName('CAPABILITY');
            var uri;
            for (var ci = 0; ci < capXMLs.length; ++ci) {
                var capXML = capXMLs[ci];

                caps.push(capXML.getAttribute('type'));

                if (capXML.getAttribute('type') == 'das1:features') {
                    var fep = capXML.getAttribute('query_uri');
                    uri = fep.substring(0, fep.length - 'features'.length);
                }
            }

            var props = {};
            var propXMLs = versionXML.getElementsByTagName('PROP');
            for (var pi = 0; pi < propXMLs.length; ++pi) {
                pusho(props, propXMLs[pi].getAttribute('name'), propXMLs[pi].getAttribute('value'));
            }

            if (uri) {
                var source = new DASSource(uri, {
                    source_uri: sourceXML.getAttribute('uri'),
                    name: sourceXML.getAttribute('title'),
                    desc: sourceXML.getAttribute('description'),
                    coords: coords,
                    props: props,
                    capabilities: caps
                });
                sources.push(source);
            }
        }

        callback(sources);
    });
};

//
// Utility functions
//

function elementValue(element, tag) {
    var children = element.getElementsByTagName(tag);
    if (children.length > 0 && children[0].firstChild) {
        var c = children[0];
        if (c.childNodes.length == 1) {
            return c.firstChild.nodeValue;
        } else {
            var s = '';
            for (var ni = 0; ni < c.childNodes.length; ++ni) {
                s += c.childNodes[ni].nodeValue;
            }
            return s;
        }
    } else {
        return null;
    }
}

function childElementOf(element) {
    if (element.hasChildNodes()) {
        var child = element.firstChild;
        do {
            if (child.nodeType == Node.ELEMENT_NODE) {
                return child;
            }
            child = child.nextSibling;
        } while (child != null);
    }
    return null;
}

function dasLinksOf(element) {
    var links = new Array();
    var maybeLinkChilden = element.getElementsByTagName('LINK');
    for (var ci = 0; ci < maybeLinkChilden.length; ++ci) {
        var linkXML = maybeLinkChilden[ci];
        if (linkXML.parentNode == element) {
            links.push(new DASLink(linkXML.firstChild ? linkXML.firstChild.nodeValue : 'Unknown', linkXML.getAttribute('href')));
        }
    }

    return links;
}

function dasNotesOf(element) {
    var notes = [];
    var maybeNotes = element.getElementsByTagName('NOTE');
    for (var ni = 0; ni < maybeNotes.length; ++ni) {
        if (maybeNotes[ni].firstChild) {
            notes.push(maybeNotes[ni].firstChild.nodeValue);
        }
    }
    return notes;
}

function doCrossDomainRequest(url, handler, credentials, custAuth) {
    // TODO: explicit error handlers?

    if (window.XDomainRequest) {
        var req = new XDomainRequest();
        req.onload = function () {
            var dom = new ActiveXObject("Microsoft.XMLDOM");
            dom.async = false;
            dom.loadXML(req.responseText);
            handler(dom);
        };
        req.open("get", url);
        req.send();
    } else {
        try {
            var req = new XMLHttpRequest();
            var timeout = setTimeout(function () {
                console.log('timing out ' + url);
                req.abort();
                handler(null, req);
            }, 5000);

            req.ontimeout = function () {
                console.log('timeout on ' + url);
            };

            req.onreadystatechange = function () {
                if (req.readyState == 4) {
                    clearTimeout(timeout);
                    if (req.status >= 200 || req.status == 0) {
                        handler(req.responseXML, req);
                    }
                }
            };
            req.open("get", url, true);
            // IE10/11 fix: The timeout property may be set only in the time interval between a call to the open method
            //              and the first call to the send method.
            req.timeout = 5000;
            if (credentials) {
                req.withCredentials = true;
            }
            if (custAuth) {
                req.setRequestHeader('X-DAS-Authorisation', custAuth);
            }
            req.overrideMimeType('text/xml');
            req.setRequestHeader('Accept', 'application/xml,*/*');
            req.send();
        } catch (e) {
            handler(null, req, e);
        }
    }
}

DASSource.prototype.doCrossDomainRequest = function (url, handler, errHandler) {
    var custAuth;
    if (this.xUser) {
        custAuth = 'Basic ' + btoa(this.xUser + ':' + this.xPass);
    }

    try {
        return doCrossDomainRequest(url, handler, this.credentials, custAuth);
    } catch (err) {
        if (errHandler) {
            errHandler(err);
        } else {
            throw err;
        }
    }
};

function isDasBooleanTrue(s) {
    s = ('' + s).toLowerCase();
    return s === 'yes' || s === 'true';
}

function isDasBooleanNotFalse(s) {
    if (!s) return false;

    s = ('' + s).toLowerCase();
    return s !== 'no' || s !== 'false';
}

function copyStylesheet(ss) {
    var nss = shallowCopy(ss);
    nss.styles = [];
    for (var si = 0; si < ss.styles.length; ++si) {
        var sh = nss.styles[si] = shallowCopy(ss.styles[si]);
        sh._methodRE = sh._labelRE = sh._typeRE = undefined;
        sh.style = shallowCopy(sh.style);
        sh.style.id = undefined;
        sh.style._gradient = undefined;
    }
    return nss;
}

if (typeof module !== 'undefined') {
    module.exports = {
        DASGroup: DASGroup,
        DASFeature: DASFeature,
        DASStylesheet: DASStylesheet,
        DASStyle: DASStyle,
        DASSource: DASSource,
        DASSegment: DASSegment,
        DASRegistry: DASRegistry,
        DASSequence: DASSequence,
        DASLink: DASLink,

        isDasBooleanTrue: isDasBooleanTrue,
        isDasBooleanNotFalse: isDasBooleanNotFalse,
        copyStylesheet: copyStylesheet,
        coordsMatch: coordsMatch
    };
}

},{"./color":4,"./utils":11}],6:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// encode.js: interface for ENCODE DCC services
//

"use strict";

var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

if (typeof require !== 'undefined') {
    var Promise = require('es6-promise').Promise;
}

function lookupEncodeURI(uri, json) {
    if (uri.indexOf('?') < 0) uri = uri + '?soft=true';

    return new Promise(function (accept, reject) {
        var req = new XMLHttpRequest();
        req.onreadystatechange = function () {
            if (req.readyState == 4) {
                if (req.status >= 300) {
                    reject('Error code ' + req.status);
                } else {
                    var resp = JSON.parse(req.response);
                    accept(json ? resp : resp.location);
                }
            }
        };

        req.open('GET', uri, true);
        req.setRequestHeader('Accept', 'application/json');
        req.responseType = 'text';
        req.send();
    });
}

function EncodeURLHolder(url) {
    this.rawurl = url;
}

EncodeURLHolder.prototype.getURLPromise = function () {
    if (this.urlPromise && this.urlPromiseValidity > Date.now()) {
        return this.urlPromise;
    } else {
        this.urlPromise = lookupEncodeURI(this.rawurl, true).then(function (resp) {
            return resp.location;
        });
        this.urlPromiseValidity = Date.now() + 12 * 3600 * 1000;
        return this.urlPromise;
    }
};

function EncodeFetchable(url, start, end, opts) {
    if (!opts) {
        if ((typeof start === 'undefined' ? 'undefined' : _typeof(start)) === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = typeof url === 'string' ? new EncodeURLHolder(url) : url;
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}

EncodeFetchable.prototype.slice = function (s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start,
        ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new EncodeFetchable(this.url, ns, ne, this.opts);
};

EncodeFetchable.prototype.fetchAsText = function (callback) {
    var self = this;
    var req = new XMLHttpRequest();
    var length;
    self.url.getURLPromise().then(function (url) {
        req.open('GET', url, true);

        if (self.end) {
            if (self.end - self.start > 100000000) {
                throw 'Monster fetch!';
            }
            req.setRequestHeader('Range', 'bytes=' + self.start + '-' + self.end);
            length = self.end - self.start + 1;
        }

        req.onreadystatechange = function () {
            if (req.readyState == 4) {
                if (req.status == 200 || req.status == 206) {
                    return callback(req.responseText);
                } else {
                    return callback(null);
                }
            }
        };
        if (self.opts.credentials) {
            req.withCredentials = true;
        }
        req.send();
    }).catch(function (err) {
        console.log(err);
        return callback(null);
    });
};

EncodeFetchable.prototype.salted = function () {
    return this;
};

EncodeFetchable.prototype.fetch = function (callback, attempt, truncatedLength) {
    var self = this;

    attempt = attempt || 1;
    if (attempt > 3) {
        return callback(null);
    }

    self.url.getURLPromise().then(function (url) {
        var req = new XMLHttpRequest();
        var length;
        req.open('GET', url, true);
        req.overrideMimeType('text/plain; charset=x-user-defined');
        if (self.end) {
            if (self.end - self.start > 100000000) {
                throw 'Monster fetch!';
            }
            req.setRequestHeader('Range', 'bytes=' + self.start + '-' + self.end);
            length = self.end - self.start + 1;
        }
        req.responseType = 'arraybuffer';
        req.onreadystatechange = function () {
            if (req.readyState == 4) {
                if (req.status == 200 || req.status == 206) {
                    if (req.response) {
                        var bl = req.response.byteLength;
                        if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                            return self.fetch(callback, attempt + 1, bl);
                        } else {
                            return callback(req.response);
                        }
                    } else if (req.mozResponseArrayBuffer) {
                        return callback(req.mozResponseArrayBuffer);
                    } else {
                        var r = req.responseText;
                        if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                            return self.fetch(callback, attempt + 1, r.length);
                        } else {
                            return callback(bstringToBuffer(req.responseText));
                        }
                    }
                } else {
                    return self.fetch(callback, attempt + 1);
                }
            }
        };
        if (self.opts.credentials) {
            req.withCredentials = true;
        }
        req.send();
    }).catch(function (err) {
        console.log(err);
    });
};

function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

if (typeof module !== 'undefined') {
    module.exports = {
        lookupEncodeURI: lookupEncodeURI,
        EncodeFetchable: EncodeFetchable
    };
}

},{"es6-promise":12}],7:[function(require,module,exports){
(function (global){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// fetchworker.js
//

"use strict";

var bin = require('./bin');
var bam = require('./bam');
var bigwig = require('./bigwig');
var encode = require('./encode');
var utils = require('./utils');

var Promise = require('es6-promise').Promise;

var connections = {};
var resolveResolvers = {};

var idSeed = 0;

global.newID = function () {
    return 'cn' + ++idSeed;
};

postMessage({ tag: 'init' });

self.onmessage = function (event) {
    var d = event.data;
    var command = event.data.command;
    var tag = event.data.tag;

    if (!command) {
        var rr = resolveResolvers[tag];
        if (rr) {
            if (d.err) {
                rr.reject(d.err);
            } else {
                rr.resolve(d.url);
            }

            delete resolveResolvers[tag];
        }
    } else if (command === 'connectBAM') {
        var id = newID();
        var resolver;
        if (d.resolver) {
            resolver = proxyResolver(d.resolver);
        }
        var bamF, baiF, indexChunks;
        if (d.blob) {
            bamF = new bin.BlobFetchable(d.blob);
            baiF = new bin.BlobFetchable(d.indexBlob);
        } else {
            bamF = new bin.URLFetchable(d.uri, { credentials: d.credentials, resolver: resolver });
            baiF = new bin.URLFetchable(d.indexUri, { credentials: d.credentials, resolver: resolver });
            indexChunks = d.indexChunks;
        }

        bam.makeBam(bamF, baiF, indexChunks, function (bamObj, err) {
            if (bamObj) {
                connections[id] = new BAMWorkerFetcher(bamObj);
                postMessage({ tag: tag, result: id });
            } else {
                postMessage({ tag: tag, error: err || "Couldn't fetch BAM" });
            }
        });
    } else if (command === 'connectBBI') {
        var id = newID();
        var resolver;
        if (d.resolver) {
            resolver = proxyResolver(d.resolver);
        }
        var bbi;
        if (d.blob) {
            bbi = new bin.BlobFetchable(d.blob);
        } else if (d.transport == 'encode') {
            bbi = new encode.EncodeFetchable(d.uri, { credentials: d.credentials });
        } else {
            bbi = new bin.URLFetchable(d.uri, { credentials: d.credentials, resolver: resolver });
        }

        bigwig.makeBwg(bbi, function (bwg, err) {
            if (bwg) {
                connections[id] = new BBIWorkerFetcher(bwg);
                postMessage({ tag: tag, result: id });
            } else {
                postMessage({ tag: tag, error: err || "Couldn't fetch BBI" });
            }
        }, d.uri);
    } else if (command === 'textxhr') {
        utils.textXHR(d.uri, function (resp, err) {
            if (resp) {
                postMessage({ tag: tag, result: resp });
            } else {
                postMessage({ tag: tag, err: err || "Couldn't fetch resource" });
            }
        });
    } else if (command === 'fetch') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({ tag: tag, error: 'No such connection: ' + event.data.connection });
        }

        con.fetch(d.tag, d.chr, d.min, d.max, d.zoom, d.opts);
    } else if (command === 'leap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({ tag: tag, error: 'No such connection: ' + event.data.connection });
        }

        con.leap(d.tag, d.chr, d.pos, d.dir);
    } else if (command === 'quantLeap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({ tag: tag, error: 'No such connection: ' + event.data.connection });
        }

        con.quantLeap(d.tag, d.chr, d.pos, d.dir, d.threshold, d.under);
    } else if (command === 'meta') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({ tag: tag, error: 'No such connection: ' + event.data.connection });
        }

        con.meta(d.tag);
    } else if (command === 'search') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({ tag: tag, error: 'No such connection: ' + event.data.connection });
        }

        con.search(d.tag, d.query, d.index);
    } else if (command === 'date') {
        return postMessage({ tag: tag, result: Date.now() | 0 });
    } else {
        postMessage({ tag: tag, error: 'Bad command ' + command });
    }
};

function BAMWorkerFetcher(bam) {
    this.bam = bam;
}

BAMWorkerFetcher.prototype.fetch = function (tag, chr, min, max, zoom, opts) {
    opts = opts || {};
    this.bam.fetch(chr, min, max, function (records, err) {
        if (records) {
            postMessage({ tag: tag, result: records, time: Date.now() | 0 });
        } else {
            postMessage({ tag: tag, error: err });
        }
    }, opts);
};

function BBIWorkerFetcher(bbi) {
    this.bbi = bbi;
}

BBIWorkerFetcher.prototype.fetch = function (tag, chr, min, max, zoom) {
    if (typeof zoom !== 'number') zoom = -1;

    var data;
    if (zoom < 0) {
        data = this.bbi.getUnzoomedView();
    } else {
        data = this.bbi.getZoomedView(zoom);
    }

    data.readWigData(chr, min, max, function (features) {
        postMessage({ tag: tag, result: features });
    });
};

BBIWorkerFetcher.prototype.meta = function (tag) {
    var scales = [1];
    for (var z = 0; z < this.bbi.zoomLevels.length; ++z) {
        scales.push(this.bbi.zoomLevels[z].reduction);
    }

    var thisB = this;
    var meta = { type: this.bbi.type,
        zoomLevels: scales,
        fieldCount: this.bbi.fieldCount,
        definedFieldCount: this.bbi.definedFieldCount,
        schema: this.bbi.schema };
    if (this.bbi.type === 'bigbed') {
        this.bbi.getExtraIndices(function (ei) {
            if (ei) {
                thisB.extraIndices = ei;
                meta.extraIndices = ei.map(function (i) {
                    return i.field;
                });
            }
            postMessage({ tag: tag, result: meta });
        });
    } else {
        postMessage({ tag: tag, result: meta });
    }
};

BBIWorkerFetcher.prototype.leap = function (tag, chr, pos, dir) {
    this.bbi.getUnzoomedView().getFirstAdjacent(chr, pos, dir, function (result, err) {
        postMessage({ tag: tag, result: result, error: err });
    });
};

BBIWorkerFetcher.prototype.quantLeap = function (tag, chr, pos, dir, threshold, under) {
    this.bbi.thresholdSearch(chr, pos, dir, threshold, function (result, err) {
        postMessage({ tag: tag, result: result, error: err });
    });
};

BBIWorkerFetcher.prototype.search = function (tag, query, index) {
    var is = this.extraIndices[0];
    is.lookup(query, function (result, err) {
        postMessage({ tag: tag, result: result, error: err });
    });
};

function proxyResolver(tag) {
    return function (url) {
        var lid = newID();
        return new Promise(function (resolve, reject) {
            resolveResolvers[lid] = { resolve: resolve, reject: reject };
            postMessage({ tag: lid,
                cmd: 'resolve',
                resolver: tag,
                url: url });
        });
    };
}

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})

},{"./bam":1,"./bigwig":2,"./bin":3,"./encode":6,"./utils":11,"es6-promise":12}],8:[function(require,module,exports){
'use strict';

/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// lh3utils.js: common support for lh3's file formats
//

if (typeof require !== 'undefined') {
    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

function Vob(b, o) {
    this.block = b;
    this.offset = o;
}

Vob.prototype.toString = function () {
    return '' + this.block + ':' + this.offset;
};

function readVob(ba, offset, allowZero) {
    var block = (ba[offset + 6] & 0xff) * 0x100000000 + (ba[offset + 5] & 0xff) * 0x1000000 + (ba[offset + 4] & 0xff) * 0x10000 + (ba[offset + 3] & 0xff) * 0x100 + (ba[offset + 2] & 0xff);
    var bint = ba[offset + 1] << 8 | ba[offset];
    if (block == 0 && bint == 0 && !allowZero) {
        return null; // Should only happen in the linear index?
    } else {
        return new Vob(block, bint);
    }
}

function unbgzf(data, lim) {
    lim = Math.min(lim || 1, data.byteLength - 50);
    var oBlockList = [];
    var ptr = [0];
    var totalSize = 0;

    while (ptr[0] < lim) {
        var ba = new Uint8Array(data, ptr[0], 12); // FIXME is this enough for all credible BGZF block headers?
        var xlen = ba[11] << 8 | ba[10];
        // dlog('xlen[' + (ptr[0]) +']=' + xlen);
        var unc = jszlib_inflate_buffer(data, 12 + xlen + ptr[0], Math.min(65536, data.byteLength - 12 - xlen - ptr[0]), ptr);
        ptr[0] += 8;
        totalSize += unc.byteLength;
        oBlockList.push(unc);
    }

    if (oBlockList.length == 1) {
        return oBlockList[0];
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = new Uint8Array(oBlockList[i]);
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

function Chunk(minv, maxv) {
    this.minv = minv;this.maxv = maxv;
}

//
// Binning (transliterated from SAM1.3 spec)
//

/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
function reg2bin(beg, end) {
    --end;
    if (beg >> 14 == end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
    if (beg >> 17 == end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
    if (beg >> 20 == end >> 20) return ((1 << 9) - 1) / 7 + (beg >> 20);
    if (beg >> 23 == end >> 23) return ((1 << 6) - 1) / 7 + (beg >> 23);
    if (beg >> 26 == end >> 26) return ((1 << 3) - 1) / 7 + (beg >> 26);
    return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
var MAX_BIN = ((1 << 18) - 1) / 7;
function reg2bins(beg, end) {
    var i = 0,
        k,
        list = [];
    --end;
    list.push(0);
    for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) {
        list.push(k);
    }for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) {
        list.push(k);
    }for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) {
        list.push(k);
    }for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) {
        list.push(k);
    }for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) {
        list.push(k);
    }return list;
}

if (typeof module !== 'undefined') {
    module.exports = {
        unbgzf: unbgzf,
        readVob: readVob,
        reg2bin: reg2bin,
        reg2bins: reg2bins,
        Chunk: Chunk
    };
}

},{"jszlib":13}],9:[function(require,module,exports){
/*
 * A JavaScript implementation of the Secure Hash Algorithm, SHA-1, as defined
 * in FIPS 180-1
 * Version 2.2 Copyright Paul Johnston 2000 - 2009.
 * Other contributors: Greg Holt, Andrew Kepert, Ydnar, Lostinet
 * Distributed under the BSD License
 * See http://pajhome.org.uk/crypt/md5 for details.
 */

"use strict";

/*
 * Configurable variables. You may need to tweak these to be compatible with
 * the server-side, but the defaults work in most cases.
 */

var hexcase = 0; /* hex output format. 0 - lowercase; 1 - uppercase        */
var b64pad = ""; /* base-64 pad character. "=" for strict RFC compliance   */

/*
 * These are the functions you'll usually want to call
 * They take string arguments and return either hex or base-64 encoded strings
 */
function hex_sha1(s) {
  return rstr2hex(rstr_sha1(str2rstr_utf8(s)));
}
function b64_sha1(s) {
  return rstr2b64(rstr_sha1(str2rstr_utf8(s)));
}
function any_sha1(s, e) {
  return rstr2any(rstr_sha1(str2rstr_utf8(s)), e);
}
function hex_hmac_sha1(k, d) {
  return rstr2hex(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)));
}
function b64_hmac_sha1(k, d) {
  return rstr2b64(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)));
}
function any_hmac_sha1(k, d, e) {
  return rstr2any(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)), e);
}

/*
 * Perform a simple self-test to see if the VM is working
 */
function sha1_vm_test() {
  return hex_sha1("abc").toLowerCase() == "a9993e364706816aba3e25717850c26c9cd0d89d";
}

/*
 * Calculate the SHA1 of a raw string
 */
function rstr_sha1(s) {
  return binb2rstr(binb_sha1(rstr2binb(s), s.length * 8));
}

/*
 * Calculate the HMAC-SHA1 of a key and some data (raw strings)
 */
function rstr_hmac_sha1(key, data) {
  var bkey = rstr2binb(key);
  if (bkey.length > 16) bkey = binb_sha1(bkey, key.length * 8);

  var ipad = Array(16),
      opad = Array(16);
  for (var i = 0; i < 16; i++) {
    ipad[i] = bkey[i] ^ 0x36363636;
    opad[i] = bkey[i] ^ 0x5C5C5C5C;
  }

  var hash = binb_sha1(ipad.concat(rstr2binb(data)), 512 + data.length * 8);
  return binb2rstr(binb_sha1(opad.concat(hash), 512 + 160));
}

/*
 * Convert a raw string to a hex string
 */
function rstr2hex(input) {
  // try { hexcase } catch(e) { hexcase=0; }
  var hex_tab = hexcase ? "0123456789ABCDEF" : "0123456789abcdef";
  var output = "";
  var x;
  for (var i = 0; i < input.length; i++) {
    x = input.charCodeAt(i);
    output += hex_tab.charAt(x >>> 4 & 0x0F) + hex_tab.charAt(x & 0x0F);
  }
  return output;
}

/*
 * Convert a raw string to a base-64 string
 */
function rstr2b64(input) {
  // try { b64pad } catch(e) { b64pad=''; }
  var tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var output = "";
  var len = input.length;
  for (var i = 0; i < len; i += 3) {
    var triplet = input.charCodeAt(i) << 16 | (i + 1 < len ? input.charCodeAt(i + 1) << 8 : 0) | (i + 2 < len ? input.charCodeAt(i + 2) : 0);
    for (var j = 0; j < 4; j++) {
      if (i * 8 + j * 6 > input.length * 8) output += b64pad;else output += tab.charAt(triplet >>> 6 * (3 - j) & 0x3F);
    }
  }
  return output;
}

/*
 * Convert a raw string to an arbitrary string encoding
 */
function rstr2any(input, encoding) {
  var divisor = encoding.length;
  var remainders = Array();
  var i, q, x, quotient;

  /* Convert to an array of 16-bit big-endian values, forming the dividend */
  var dividend = Array(Math.ceil(input.length / 2));
  for (i = 0; i < dividend.length; i++) {
    dividend[i] = input.charCodeAt(i * 2) << 8 | input.charCodeAt(i * 2 + 1);
  }

  /*
   * Repeatedly perform a long division. The binary array forms the dividend,
   * the length of the encoding is the divisor. Once computed, the quotient
   * forms the dividend for the next step. We stop when the dividend is zero.
   * All remainders are stored for later use.
   */
  while (dividend.length > 0) {
    quotient = Array();
    x = 0;
    for (i = 0; i < dividend.length; i++) {
      x = (x << 16) + dividend[i];
      q = Math.floor(x / divisor);
      x -= q * divisor;
      if (quotient.length > 0 || q > 0) quotient[quotient.length] = q;
    }
    remainders[remainders.length] = x;
    dividend = quotient;
  }

  /* Convert the remainders to the output string */
  var output = "";
  for (i = remainders.length - 1; i >= 0; i--) {
    output += encoding.charAt(remainders[i]);
  } /* Append leading zero equivalents */
  var full_length = Math.ceil(input.length * 8 / (Math.log(encoding.length) / Math.log(2)));
  for (i = output.length; i < full_length; i++) {
    output = encoding[0] + output;
  }return output;
}

/*
 * Encode a string as utf-8.
 * For efficiency, this assumes the input is valid utf-16.
 */
function str2rstr_utf8(input) {
  var output = "";
  var i = -1;
  var x, y;

  while (++i < input.length) {
    /* Decode utf-16 surrogate pairs */
    x = input.charCodeAt(i);
    y = i + 1 < input.length ? input.charCodeAt(i + 1) : 0;
    if (0xD800 <= x && x <= 0xDBFF && 0xDC00 <= y && y <= 0xDFFF) {
      x = 0x10000 + ((x & 0x03FF) << 10) + (y & 0x03FF);
      i++;
    }

    /* Encode output as utf-8 */
    if (x <= 0x7F) output += String.fromCharCode(x);else if (x <= 0x7FF) output += String.fromCharCode(0xC0 | x >>> 6 & 0x1F, 0x80 | x & 0x3F);else if (x <= 0xFFFF) output += String.fromCharCode(0xE0 | x >>> 12 & 0x0F, 0x80 | x >>> 6 & 0x3F, 0x80 | x & 0x3F);else if (x <= 0x1FFFFF) output += String.fromCharCode(0xF0 | x >>> 18 & 0x07, 0x80 | x >>> 12 & 0x3F, 0x80 | x >>> 6 & 0x3F, 0x80 | x & 0x3F);
  }
  return output;
}

/*
 * Encode a string as utf-16
 */
function str2rstr_utf16le(input) {
  var output = "";
  for (var i = 0; i < input.length; i++) {
    output += String.fromCharCode(input.charCodeAt(i) & 0xFF, input.charCodeAt(i) >>> 8 & 0xFF);
  }return output;
}

function str2rstr_utf16be(input) {
  var output = "";
  for (var i = 0; i < input.length; i++) {
    output += String.fromCharCode(input.charCodeAt(i) >>> 8 & 0xFF, input.charCodeAt(i) & 0xFF);
  }return output;
}

/*
 * Convert a raw string to an array of big-endian words
 * Characters >255 have their high-byte silently ignored.
 */
function rstr2binb(input) {
  var output = Array(input.length >> 2);
  for (var i = 0; i < output.length; i++) {
    output[i] = 0;
  }for (var i = 0; i < input.length * 8; i += 8) {
    output[i >> 5] |= (input.charCodeAt(i / 8) & 0xFF) << 24 - i % 32;
  }return output;
}

/*
 * Convert an array of big-endian words to a string
 */
function binb2rstr(input) {
  var output = "";
  for (var i = 0; i < input.length * 32; i += 8) {
    output += String.fromCharCode(input[i >> 5] >>> 24 - i % 32 & 0xFF);
  }return output;
}

/*
 * Calculate the SHA-1 of an array of big-endian words, and a bit length
 */
function binb_sha1(x, len) {
  /* append padding */
  x[len >> 5] |= 0x80 << 24 - len % 32;
  x[(len + 64 >> 9 << 4) + 15] = len;

  var w = Array(80);
  var a = 1732584193;
  var b = -271733879;
  var c = -1732584194;
  var d = 271733878;
  var e = -1009589776;

  for (var i = 0; i < x.length; i += 16) {
    var olda = a;
    var oldb = b;
    var oldc = c;
    var oldd = d;
    var olde = e;

    for (var j = 0; j < 80; j++) {
      if (j < 16) w[j] = x[i + j];else w[j] = bit_rol(w[j - 3] ^ w[j - 8] ^ w[j - 14] ^ w[j - 16], 1);
      var t = safe_add(safe_add(bit_rol(a, 5), sha1_ft(j, b, c, d)), safe_add(safe_add(e, w[j]), sha1_kt(j)));
      e = d;
      d = c;
      c = bit_rol(b, 30);
      b = a;
      a = t;
    }

    a = safe_add(a, olda);
    b = safe_add(b, oldb);
    c = safe_add(c, oldc);
    d = safe_add(d, oldd);
    e = safe_add(e, olde);
  }
  return Array(a, b, c, d, e);
}

/*
 * Perform the appropriate triplet combination function for the current
 * iteration
 */
function sha1_ft(t, b, c, d) {
  if (t < 20) return b & c | ~b & d;
  if (t < 40) return b ^ c ^ d;
  if (t < 60) return b & c | b & d | c & d;
  return b ^ c ^ d;
}

/*
 * Determine the appropriate additive constant for the current iteration
 */
function sha1_kt(t) {
  return t < 20 ? 1518500249 : t < 40 ? 1859775393 : t < 60 ? -1894007588 : -899497514;
}

/*
 * Add integers, wrapping at 2^32. This uses 16-bit operations internally
 * to work around bugs in some JS interpreters.
 */
function safe_add(x, y) {
  var lsw = (x & 0xFFFF) + (y & 0xFFFF);
  var msw = (x >> 16) + (y >> 16) + (lsw >> 16);
  return msw << 16 | lsw & 0xFFFF;
}

/*
 * Bitwise rotate a 32-bit number to the left.
 */
function bit_rol(num, cnt) {
  return num << cnt | num >>> 32 - cnt;
}

if (typeof module !== 'undefined') {
  module.exports = {
    b64_sha1: b64_sha1,
    hex_sha1: hex_sha1
  };
}

},{}],10:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// spans.js: JavaScript Intset/Location port.
//

"use strict";

function Range(min, max) {
    if (typeof min != 'number' || typeof max != 'number') throw 'Bad range ' + min + ',' + max;
    this._min = min;
    this._max = max;
}

Range.prototype.min = function () {
    return this._min;
};

Range.prototype.max = function () {
    return this._max;
};

Range.prototype.contains = function (pos) {
    return pos >= this._min && pos <= this._max;
};

Range.prototype.isContiguous = function () {
    return true;
};

Range.prototype.ranges = function () {
    return [this];
};

Range.prototype._pushRanges = function (ranges) {
    ranges.push(this);
};

Range.prototype.toString = function () {
    return '[' + this._min + '-' + this._max + ']';
};

function _Compound(ranges) {
    // given: a set of unsorted possibly overlapping ranges
    // sort the input ranges
    var sorted = ranges.sort(_rangeOrder);
    // merge overlaps between adjacent ranges
    var merged = [];
    var current = sorted.shift();
    sorted.forEach(function (range) {
        if (range._min <= current._max) {
            if (range._max > current._max) {
                current._max = range._max;
            }
        } else {
            merged.push(current);
            current = range;
        }
    });
    merged.push(current);
    this._ranges = merged;
}

_Compound.prototype.min = function () {
    return this._ranges[0].min();
};

_Compound.prototype.max = function () {
    return this._ranges[this._ranges.length - 1].max();
};

// returns the index of the first range that is not less than pos
_Compound.prototype.lower_bound = function (pos) {
    // first check if pos is out of range
    var r = this.ranges();
    if (pos > this.max()) return r.length;
    if (pos < this.min()) return 0;
    // do a binary search
    var a = 0,
        b = r.length - 1;
    while (a <= b) {
        var m = Math.floor((a + b) / 2);
        if (pos > r[m]._max) {
            a = m + 1;
        } else if (pos < r[m]._min) {
            b = m - 1;
        } else {
            return m;
        }
    }
    return a;
};

_Compound.prototype.contains = function (pos) {
    var lb = this.lower_bound(pos);
    if (lb < this._ranges.length && this._ranges[lb].contains(pos)) {
        return true;
    }
    return false;
};

_Compound.prototype.insertRange = function (range) {
    var lb = this.lower_bound(range._min);
    if (lb === this._ranges.length) {
        // range follows this
        this._ranges.push(range);
        return;
    }

    var r = this.ranges();
    if (range._max < r[lb]._min) {
        // range preceeds lb
        this._ranges.splice(lb, 0, range);
        return;
    }

    // range overlaps lb (at least)
    if (r[lb]._min < range._min) range._min = r[lb]._min;
    var ub = lb + 1;
    while (ub < r.length && r[ub]._min <= range._max) {
        ub++;
    }
    ub--;
    // ub is the upper bound of the new range
    if (r[ub]._max > range._max) range._max = r[ub]._max;

    // splice range into this._ranges
    this._ranges.splice(lb, ub - lb + 1, range);
    return;
};

_Compound.prototype.isContiguous = function () {
    return this._ranges.length > 1;
};

_Compound.prototype.ranges = function () {
    return this._ranges;
};

_Compound.prototype._pushRanges = function (ranges) {
    for (var ri = 0; ri < this._ranges.length; ++ri) {
        ranges.push(this._ranges[ri]);
    }
};

_Compound.prototype.toString = function () {
    var s = '';
    for (var r = 0; r < this._ranges.length; ++r) {
        if (r > 0) {
            s = s + ',';
        }
        s = s + this._ranges[r].toString();
    }
    return s;
};

function union(s0, s1) {
    if (!(s0 instanceof _Compound)) {
        if (!(s0 instanceof Array)) s0 = [s0];
        s0 = new _Compound(s0);
    }

    if (s1) s0.insertRange(s1);

    return s0;
}

function intersection(s0, s1) {
    var r0 = s0.ranges();
    var r1 = s1.ranges();
    var l0 = r0.length,
        l1 = r1.length;
    var i0 = 0,
        i1 = 0;
    var or = [];

    while (i0 < l0 && i1 < l1) {
        var s0 = r0[i0],
            s1 = r1[i1];
        var lapMin = Math.max(s0.min(), s1.min());
        var lapMax = Math.min(s0.max(), s1.max());
        if (lapMax >= lapMin) {
            or.push(new Range(lapMin, lapMax));
        }
        if (s0.max() > s1.max()) {
            ++i1;
        } else {
            ++i0;
        }
    }

    if (or.length == 0) {
        return null; // FIXME
    } else if (or.length == 1) {
        return or[0];
    } else {
        return new _Compound(or);
    }
}

function coverage(s) {
    var tot = 0;
    var rl = s.ranges();
    for (var ri = 0; ri < rl.length; ++ri) {
        var r = rl[ri];
        tot += r.max() - r.min() + 1;
    }
    return tot;
}

function rangeOrder(a, b) {
    if (a.min() < b.min()) {
        return -1;
    } else if (a.min() > b.min()) {
        return 1;
    } else if (a.max() < b.max()) {
        return -1;
    } else if (b.max() > a.max()) {
        return 1;
    } else {
        return 0;
    }
}

function _rangeOrder(a, b) {
    if (a._min < b._min) {
        return -1;
    } else if (a._min > b._min) {
        return 1;
    } else if (a._max < b._max) {
        return -1;
    } else if (b._max > a._max) {
        return 1;
    } else {
        return 0;
    }
}

if (typeof module !== 'undefined') {
    module.exports = {
        Range: Range,
        union: union,
        intersection: intersection,
        coverage: coverage,
        rangeOver: rangeOrder,
        _rangeOrder: _rangeOrder
    };
}

},{}],11:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// utils.js: odds, sods, and ends.
//

"use strict";

var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

if (typeof require !== 'undefined') {
    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;
}

var NUM_REGEXP = new RegExp('[0-9]+');

function stringToNumbersArray(str) {
    var nums = new Array();
    var m;
    while (m = NUM_REGEXP.exec(str)) {
        nums.push(m[0]);
        str = str.substring(m.index + m[0].length);
    }
    return nums;
}

var STRICT_NUM_REGEXP = new RegExp('^[0-9]+$');

function stringToInt(str) {
    str = str.replace(new RegExp(',', 'g'), '');
    if (!STRICT_NUM_REGEXP.test(str)) {
        return null;
    }
    return str | 0;
}

function pushnew(a, v) {
    for (var i = 0; i < a.length; ++i) {
        if (a[i] == v) {
            return;
        }
    }
    a.push(v);
}

function pusho(obj, k, v) {
    if (obj[k]) {
        obj[k].push(v);
    } else {
        obj[k] = [v];
    }
}

function pushnewo(obj, k, v) {
    var a = obj[k];
    if (a) {
        for (var i = 0; i < a.length; ++i) {
            // indexOf requires JS16 :-(.
            if (a[i] == v) {
                return;
            }
        }
        a.push(v);
    } else {
        obj[k] = [v];
    }
}

function pick(a, b, c, d) {
    if (a) {
        return a;
    } else if (b) {
        return b;
    } else if (c) {
        return c;
    } else if (d) {
        return d;
    }
}

function pushnew(l, o) {
    for (var i = 0; i < l.length; ++i) {
        if (l[i] == o) {
            return;
        }
    }
    l.push(o);
}

function arrayEquals(a, b) {
    if (!a || !b || a.length !== b.length) return false;

    for (var i = 0; i < a.length; ++i) {
        if (a[i] !== b[i]) return false;
    }

    return true;
}

function arrayIndexOf(a, x) {
    if (!a) {
        return -1;
    }

    for (var i = 0; i < a.length; ++i) {
        if (a[i] === x) {
            return i;
        }
    }
    return -1;
}

function arrayRemove(a, x) {
    var i = arrayIndexOf(a, x);
    if (i >= 0) {
        a.splice(i, 1);
        return true;
    }
    return false;
}

//
// DOM utilities
//


function makeElement(tag, children, attribs, styles) {
    var ele = document.createElement(tag);
    if (children) {
        if (!(children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (c) {
                if (typeof c == 'string') {
                    c = document.createTextNode(c);
                } else if (typeof c == 'number') {
                    c = document.createTextNode('' + c);
                }
                ele.appendChild(c);
            }
        }
    }

    if (attribs) {
        for (var l in attribs) {
            try {
                ele[l] = attribs[l];
            } catch (e) {
                console.log('error setting ' + l);
                throw e;
            }
        }
    }
    if (styles) {
        for (var l in styles) {
            ele.style[l] = styles[l];
        }
    }
    return ele;
}

function makeElementNS(namespace, tag, children, attribs) {
    var ele = document.createElementNS(namespace, tag);
    if (children) {
        if (!(children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (typeof c == 'string') {
                c = document.createTextNode(c);
            }
            ele.appendChild(c);
        }
    }

    setAttrs(ele, attribs);
    return ele;
}

var attr_name_cache = {};

function setAttr(node, key, value) {
    var attr = attr_name_cache[key];
    if (!attr) {
        var _attr = '';
        for (var c = 0; c < key.length; ++c) {
            var cc = key.substring(c, c + 1);
            var lcc = cc.toLowerCase();
            if (lcc != cc) {
                _attr = _attr + '-' + lcc;
            } else {
                _attr = _attr + cc;
            }
        }
        attr_name_cache[key] = _attr;
        attr = _attr;
    }
    node.setAttribute(attr, value);
}

function setAttrs(node, attribs) {
    if (attribs) {
        for (var l in attribs) {
            setAttr(node, l, attribs[l]);
        }
    }
}

function removeChildren(node) {
    if (!node || !node.childNodes) {
        return;
    }

    while (node.childNodes.length > 0) {
        node.removeChild(node.firstChild);
    }
}

//
// WARNING: not for general use!
//

function miniJSONify(o, exc) {
    if (typeof o === 'undefined') {
        return 'undefined';
    } else if (o == null) {
        return 'null';
    } else if (typeof o == 'string') {
        return "'" + o + "'";
    } else if (typeof o == 'number') {
        return "" + o;
    } else if (typeof o == 'boolean') {
        return "" + o;
    } else if ((typeof o === 'undefined' ? 'undefined' : _typeof(o)) == 'object') {
        if (o instanceof Array) {
            var s = null;
            for (var i = 0; i < o.length; ++i) {
                s = (s == null ? '' : s + ', ') + miniJSONify(o[i], exc);
            }
            return '[' + (s ? s : '') + ']';
        } else {
            exc = exc || {};
            var s = null;
            for (var k in o) {
                if (exc[k]) continue;
                if (k != undefined && typeof o[k] != 'function') {
                    s = (s == null ? '' : s + ', ') + k + ': ' + miniJSONify(o[k], exc);
                }
            }
            return '{' + (s ? s : '') + '}';
        }
    } else {
        return typeof o === 'undefined' ? 'undefined' : _typeof(o);
    }
}

function shallowCopy(o) {
    var n = {};
    for (var k in o) {
        n[k] = o[k];
    }
    return n;
}

function Observed(x) {
    this.value = x;
    this.listeners = [];
}

Observed.prototype.addListener = function (f) {
    this.listeners.push(f);
};

Observed.prototype.addListenerAndFire = function (f) {
    this.listeners.push(f);
    f(this.value);
};

Observed.prototype.removeListener = function (f) {
    arrayRemove(this.listeners, f);
};

Observed.prototype.get = function () {
    return this.value;
};

Observed.prototype.set = function (x) {
    this.value = x;
    for (var i = 0; i < this.listeners.length; ++i) {
        this.listeners[i](x);
    }
};

function Awaited() {
    this.queue = [];
}

Awaited.prototype.provide = function (x) {
    if (this.res !== undefined) {
        throw "Resource has already been provided.";
    }

    this.res = x;
    for (var i = 0; i < this.queue.length; ++i) {
        this.queue[i](x);
    }
    this.queue = null; // avoid leaking closures.
};

Awaited.prototype.await = function (f) {
    if (this.res !== undefined) {
        f(this.res);
        return this.res;
    } else {
        this.queue.push(f);
    }
};

var __dalliance_saltSeed = 0;

function saltURL(url) {
    return url + '?salt=' + b64_sha1('' + Date.now() + ',' + ++__dalliance_saltSeed);
}

function textXHR(url, callback, opts) {
    if (opts && opts.salt) url = saltURL(url);

    try {
        var timeout;
        if (opts && opts.timeout) {
            timeout = setTimeout(function () {
                console.log('timing out ' + url);
                req.abort();
                return callback(null, 'Timeout');
            }, opts.timeout);
        }

        var req = new XMLHttpRequest();
        req.onreadystatechange = function () {
            if (req.readyState == 4) {
                if (timeout) clearTimeout(timeout);
                if (req.status < 200 || req.status >= 300) {
                    callback(null, 'Error code ' + req.status);
                } else {
                    callback(req.responseText);
                }
            }
        };

        req.open('GET', url, true);
        req.responseType = 'text';

        if (opts && opts.credentials) {
            req.withCredentials = true;
        }
        req.send();
    } catch (e) {
        callback(null, 'Exception ' + e);
    }
}

function relativeURL(base, rel) {
    // FIXME quite naive -- good enough for trackhubs?

    if (rel.indexOf('http:') === 0 || rel.indexOf('https:') === 0 || rel.indexOf('ftp:') === 0) {
        return rel;
    }

    var li = base.lastIndexOf('/');
    if (li >= 0) {
        return base.substr(0, li + 1) + rel;
    } else {
        return rel;
    }
}

var AMINO_ACID_TRANSLATION = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',
    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',
    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'TAT': 'Y',
    'TAC': 'Y',
    'TAA': '*', // stop
    'TAG': '*', // stop
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'TGT': 'C',
    'TGC': 'C',
    'TGA': '*', // stop
    'TGG': 'W',
    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G'
};

function resolveUrlToPage(rel) {
    return makeElement('a', null, { href: rel }).href;
}

//
// Missing APIs
// 

if (!('trim' in String.prototype)) {
    String.prototype.trim = function () {
        return this.replace(/^\s+/, '').replace(/\s+$/, '');
    };
}

if (typeof module !== 'undefined') {
    module.exports = {
        textXHR: textXHR,
        relativeURL: relativeURL,
        resolveUrlToPage: resolveUrlToPage,
        shallowCopy: shallowCopy,
        pusho: pusho,
        pushnew: pushnew,
        pushnewo: pushnewo,
        arrayIndexOf: arrayIndexOf,
        arrayEquals: arrayEquals,
        pick: pick,

        makeElement: makeElement,
        makeElementNS: makeElementNS,
        removeChildren: removeChildren,

        miniJSONify: miniJSONify,

        Observed: Observed,
        Awaited: Awaited,

        AMINO_ACID_TRANSLATION: AMINO_ACID_TRANSLATION
    };
}

},{"./sha1":9}],12:[function(require,module,exports){
(function (process,global){
/*!
 * @overview es6-promise - a tiny implementation of Promises/A+.
 * @copyright Copyright (c) 2014 Yehuda Katz, Tom Dale, Stefan Penner and contributors (Conversion to ES6 API by Jake Archibald)
 * @license   Licensed under MIT license
 *            See https://raw.githubusercontent.com/jakearchibald/es6-promise/master/LICENSE
 * @version   3.0.2
 */

(function() {
    "use strict";
    function lib$es6$promise$utils$$objectOrFunction(x) {
      return typeof x === 'function' || (typeof x === 'object' && x !== null);
    }

    function lib$es6$promise$utils$$isFunction(x) {
      return typeof x === 'function';
    }

    function lib$es6$promise$utils$$isMaybeThenable(x) {
      return typeof x === 'object' && x !== null;
    }

    var lib$es6$promise$utils$$_isArray;
    if (!Array.isArray) {
      lib$es6$promise$utils$$_isArray = function (x) {
        return Object.prototype.toString.call(x) === '[object Array]';
      };
    } else {
      lib$es6$promise$utils$$_isArray = Array.isArray;
    }

    var lib$es6$promise$utils$$isArray = lib$es6$promise$utils$$_isArray;
    var lib$es6$promise$asap$$len = 0;
    var lib$es6$promise$asap$$toString = {}.toString;
    var lib$es6$promise$asap$$vertxNext;
    var lib$es6$promise$asap$$customSchedulerFn;

    var lib$es6$promise$asap$$asap = function asap(callback, arg) {
      lib$es6$promise$asap$$queue[lib$es6$promise$asap$$len] = callback;
      lib$es6$promise$asap$$queue[lib$es6$promise$asap$$len + 1] = arg;
      lib$es6$promise$asap$$len += 2;
      if (lib$es6$promise$asap$$len === 2) {
        // If len is 2, that means that we need to schedule an async flush.
        // If additional callbacks are queued before the queue is flushed, they
        // will be processed by this flush that we are scheduling.
        if (lib$es6$promise$asap$$customSchedulerFn) {
          lib$es6$promise$asap$$customSchedulerFn(lib$es6$promise$asap$$flush);
        } else {
          lib$es6$promise$asap$$scheduleFlush();
        }
      }
    }

    function lib$es6$promise$asap$$setScheduler(scheduleFn) {
      lib$es6$promise$asap$$customSchedulerFn = scheduleFn;
    }

    function lib$es6$promise$asap$$setAsap(asapFn) {
      lib$es6$promise$asap$$asap = asapFn;
    }

    var lib$es6$promise$asap$$browserWindow = (typeof window !== 'undefined') ? window : undefined;
    var lib$es6$promise$asap$$browserGlobal = lib$es6$promise$asap$$browserWindow || {};
    var lib$es6$promise$asap$$BrowserMutationObserver = lib$es6$promise$asap$$browserGlobal.MutationObserver || lib$es6$promise$asap$$browserGlobal.WebKitMutationObserver;
    var lib$es6$promise$asap$$isNode = typeof process !== 'undefined' && {}.toString.call(process) === '[object process]';

    // test for web worker but not in IE10
    var lib$es6$promise$asap$$isWorker = typeof Uint8ClampedArray !== 'undefined' &&
      typeof importScripts !== 'undefined' &&
      typeof MessageChannel !== 'undefined';

    // node
    function lib$es6$promise$asap$$useNextTick() {
      // node version 0.10.x displays a deprecation warning when nextTick is used recursively
      // see https://github.com/cujojs/when/issues/410 for details
      return function() {
        process.nextTick(lib$es6$promise$asap$$flush);
      };
    }

    // vertx
    function lib$es6$promise$asap$$useVertxTimer() {
      return function() {
        lib$es6$promise$asap$$vertxNext(lib$es6$promise$asap$$flush);
      };
    }

    function lib$es6$promise$asap$$useMutationObserver() {
      var iterations = 0;
      var observer = new lib$es6$promise$asap$$BrowserMutationObserver(lib$es6$promise$asap$$flush);
      var node = document.createTextNode('');
      observer.observe(node, { characterData: true });

      return function() {
        node.data = (iterations = ++iterations % 2);
      };
    }

    // web worker
    function lib$es6$promise$asap$$useMessageChannel() {
      var channel = new MessageChannel();
      channel.port1.onmessage = lib$es6$promise$asap$$flush;
      return function () {
        channel.port2.postMessage(0);
      };
    }

    function lib$es6$promise$asap$$useSetTimeout() {
      return function() {
        setTimeout(lib$es6$promise$asap$$flush, 1);
      };
    }

    var lib$es6$promise$asap$$queue = new Array(1000);
    function lib$es6$promise$asap$$flush() {
      for (var i = 0; i < lib$es6$promise$asap$$len; i+=2) {
        var callback = lib$es6$promise$asap$$queue[i];
        var arg = lib$es6$promise$asap$$queue[i+1];

        callback(arg);

        lib$es6$promise$asap$$queue[i] = undefined;
        lib$es6$promise$asap$$queue[i+1] = undefined;
      }

      lib$es6$promise$asap$$len = 0;
    }

    function lib$es6$promise$asap$$attemptVertx() {
      try {
        var r = require;
        var vertx = r('vertx');
        lib$es6$promise$asap$$vertxNext = vertx.runOnLoop || vertx.runOnContext;
        return lib$es6$promise$asap$$useVertxTimer();
      } catch(e) {
        return lib$es6$promise$asap$$useSetTimeout();
      }
    }

    var lib$es6$promise$asap$$scheduleFlush;
    // Decide what async method to use to triggering processing of queued callbacks:
    if (lib$es6$promise$asap$$isNode) {
      lib$es6$promise$asap$$scheduleFlush = lib$es6$promise$asap$$useNextTick();
    } else if (lib$es6$promise$asap$$BrowserMutationObserver) {
      lib$es6$promise$asap$$scheduleFlush = lib$es6$promise$asap$$useMutationObserver();
    } else if (lib$es6$promise$asap$$isWorker) {
      lib$es6$promise$asap$$scheduleFlush = lib$es6$promise$asap$$useMessageChannel();
    } else if (lib$es6$promise$asap$$browserWindow === undefined && typeof require === 'function') {
      lib$es6$promise$asap$$scheduleFlush = lib$es6$promise$asap$$attemptVertx();
    } else {
      lib$es6$promise$asap$$scheduleFlush = lib$es6$promise$asap$$useSetTimeout();
    }

    function lib$es6$promise$$internal$$noop() {}

    var lib$es6$promise$$internal$$PENDING   = void 0;
    var lib$es6$promise$$internal$$FULFILLED = 1;
    var lib$es6$promise$$internal$$REJECTED  = 2;

    var lib$es6$promise$$internal$$GET_THEN_ERROR = new lib$es6$promise$$internal$$ErrorObject();

    function lib$es6$promise$$internal$$selfFulfillment() {
      return new TypeError("You cannot resolve a promise with itself");
    }

    function lib$es6$promise$$internal$$cannotReturnOwn() {
      return new TypeError('A promises callback cannot return that same promise.');
    }

    function lib$es6$promise$$internal$$getThen(promise) {
      try {
        return promise.then;
      } catch(error) {
        lib$es6$promise$$internal$$GET_THEN_ERROR.error = error;
        return lib$es6$promise$$internal$$GET_THEN_ERROR;
      }
    }

    function lib$es6$promise$$internal$$tryThen(then, value, fulfillmentHandler, rejectionHandler) {
      try {
        then.call(value, fulfillmentHandler, rejectionHandler);
      } catch(e) {
        return e;
      }
    }

    function lib$es6$promise$$internal$$handleForeignThenable(promise, thenable, then) {
       lib$es6$promise$asap$$asap(function(promise) {
        var sealed = false;
        var error = lib$es6$promise$$internal$$tryThen(then, thenable, function(value) {
          if (sealed) { return; }
          sealed = true;
          if (thenable !== value) {
            lib$es6$promise$$internal$$resolve(promise, value);
          } else {
            lib$es6$promise$$internal$$fulfill(promise, value);
          }
        }, function(reason) {
          if (sealed) { return; }
          sealed = true;

          lib$es6$promise$$internal$$reject(promise, reason);
        }, 'Settle: ' + (promise._label || ' unknown promise'));

        if (!sealed && error) {
          sealed = true;
          lib$es6$promise$$internal$$reject(promise, error);
        }
      }, promise);
    }

    function lib$es6$promise$$internal$$handleOwnThenable(promise, thenable) {
      if (thenable._state === lib$es6$promise$$internal$$FULFILLED) {
        lib$es6$promise$$internal$$fulfill(promise, thenable._result);
      } else if (thenable._state === lib$es6$promise$$internal$$REJECTED) {
        lib$es6$promise$$internal$$reject(promise, thenable._result);
      } else {
        lib$es6$promise$$internal$$subscribe(thenable, undefined, function(value) {
          lib$es6$promise$$internal$$resolve(promise, value);
        }, function(reason) {
          lib$es6$promise$$internal$$reject(promise, reason);
        });
      }
    }

    function lib$es6$promise$$internal$$handleMaybeThenable(promise, maybeThenable) {
      if (maybeThenable.constructor === promise.constructor) {
        lib$es6$promise$$internal$$handleOwnThenable(promise, maybeThenable);
      } else {
        var then = lib$es6$promise$$internal$$getThen(maybeThenable);

        if (then === lib$es6$promise$$internal$$GET_THEN_ERROR) {
          lib$es6$promise$$internal$$reject(promise, lib$es6$promise$$internal$$GET_THEN_ERROR.error);
        } else if (then === undefined) {
          lib$es6$promise$$internal$$fulfill(promise, maybeThenable);
        } else if (lib$es6$promise$utils$$isFunction(then)) {
          lib$es6$promise$$internal$$handleForeignThenable(promise, maybeThenable, then);
        } else {
          lib$es6$promise$$internal$$fulfill(promise, maybeThenable);
        }
      }
    }

    function lib$es6$promise$$internal$$resolve(promise, value) {
      if (promise === value) {
        lib$es6$promise$$internal$$reject(promise, lib$es6$promise$$internal$$selfFulfillment());
      } else if (lib$es6$promise$utils$$objectOrFunction(value)) {
        lib$es6$promise$$internal$$handleMaybeThenable(promise, value);
      } else {
        lib$es6$promise$$internal$$fulfill(promise, value);
      }
    }

    function lib$es6$promise$$internal$$publishRejection(promise) {
      if (promise._onerror) {
        promise._onerror(promise._result);
      }

      lib$es6$promise$$internal$$publish(promise);
    }

    function lib$es6$promise$$internal$$fulfill(promise, value) {
      if (promise._state !== lib$es6$promise$$internal$$PENDING) { return; }

      promise._result = value;
      promise._state = lib$es6$promise$$internal$$FULFILLED;

      if (promise._subscribers.length !== 0) {
        lib$es6$promise$asap$$asap(lib$es6$promise$$internal$$publish, promise);
      }
    }

    function lib$es6$promise$$internal$$reject(promise, reason) {
      if (promise._state !== lib$es6$promise$$internal$$PENDING) { return; }
      promise._state = lib$es6$promise$$internal$$REJECTED;
      promise._result = reason;

      lib$es6$promise$asap$$asap(lib$es6$promise$$internal$$publishRejection, promise);
    }

    function lib$es6$promise$$internal$$subscribe(parent, child, onFulfillment, onRejection) {
      var subscribers = parent._subscribers;
      var length = subscribers.length;

      parent._onerror = null;

      subscribers[length] = child;
      subscribers[length + lib$es6$promise$$internal$$FULFILLED] = onFulfillment;
      subscribers[length + lib$es6$promise$$internal$$REJECTED]  = onRejection;

      if (length === 0 && parent._state) {
        lib$es6$promise$asap$$asap(lib$es6$promise$$internal$$publish, parent);
      }
    }

    function lib$es6$promise$$internal$$publish(promise) {
      var subscribers = promise._subscribers;
      var settled = promise._state;

      if (subscribers.length === 0) { return; }

      var child, callback, detail = promise._result;

      for (var i = 0; i < subscribers.length; i += 3) {
        child = subscribers[i];
        callback = subscribers[i + settled];

        if (child) {
          lib$es6$promise$$internal$$invokeCallback(settled, child, callback, detail);
        } else {
          callback(detail);
        }
      }

      promise._subscribers.length = 0;
    }

    function lib$es6$promise$$internal$$ErrorObject() {
      this.error = null;
    }

    var lib$es6$promise$$internal$$TRY_CATCH_ERROR = new lib$es6$promise$$internal$$ErrorObject();

    function lib$es6$promise$$internal$$tryCatch(callback, detail) {
      try {
        return callback(detail);
      } catch(e) {
        lib$es6$promise$$internal$$TRY_CATCH_ERROR.error = e;
        return lib$es6$promise$$internal$$TRY_CATCH_ERROR;
      }
    }

    function lib$es6$promise$$internal$$invokeCallback(settled, promise, callback, detail) {
      var hasCallback = lib$es6$promise$utils$$isFunction(callback),
          value, error, succeeded, failed;

      if (hasCallback) {
        value = lib$es6$promise$$internal$$tryCatch(callback, detail);

        if (value === lib$es6$promise$$internal$$TRY_CATCH_ERROR) {
          failed = true;
          error = value.error;
          value = null;
        } else {
          succeeded = true;
        }

        if (promise === value) {
          lib$es6$promise$$internal$$reject(promise, lib$es6$promise$$internal$$cannotReturnOwn());
          return;
        }

      } else {
        value = detail;
        succeeded = true;
      }

      if (promise._state !== lib$es6$promise$$internal$$PENDING) {
        // noop
      } else if (hasCallback && succeeded) {
        lib$es6$promise$$internal$$resolve(promise, value);
      } else if (failed) {
        lib$es6$promise$$internal$$reject(promise, error);
      } else if (settled === lib$es6$promise$$internal$$FULFILLED) {
        lib$es6$promise$$internal$$fulfill(promise, value);
      } else if (settled === lib$es6$promise$$internal$$REJECTED) {
        lib$es6$promise$$internal$$reject(promise, value);
      }
    }

    function lib$es6$promise$$internal$$initializePromise(promise, resolver) {
      try {
        resolver(function resolvePromise(value){
          lib$es6$promise$$internal$$resolve(promise, value);
        }, function rejectPromise(reason) {
          lib$es6$promise$$internal$$reject(promise, reason);
        });
      } catch(e) {
        lib$es6$promise$$internal$$reject(promise, e);
      }
    }

    function lib$es6$promise$enumerator$$Enumerator(Constructor, input) {
      var enumerator = this;

      enumerator._instanceConstructor = Constructor;
      enumerator.promise = new Constructor(lib$es6$promise$$internal$$noop);

      if (enumerator._validateInput(input)) {
        enumerator._input     = input;
        enumerator.length     = input.length;
        enumerator._remaining = input.length;

        enumerator._init();

        if (enumerator.length === 0) {
          lib$es6$promise$$internal$$fulfill(enumerator.promise, enumerator._result);
        } else {
          enumerator.length = enumerator.length || 0;
          enumerator._enumerate();
          if (enumerator._remaining === 0) {
            lib$es6$promise$$internal$$fulfill(enumerator.promise, enumerator._result);
          }
        }
      } else {
        lib$es6$promise$$internal$$reject(enumerator.promise, enumerator._validationError());
      }
    }

    lib$es6$promise$enumerator$$Enumerator.prototype._validateInput = function(input) {
      return lib$es6$promise$utils$$isArray(input);
    };

    lib$es6$promise$enumerator$$Enumerator.prototype._validationError = function() {
      return new Error('Array Methods must be provided an Array');
    };

    lib$es6$promise$enumerator$$Enumerator.prototype._init = function() {
      this._result = new Array(this.length);
    };

    var lib$es6$promise$enumerator$$default = lib$es6$promise$enumerator$$Enumerator;

    lib$es6$promise$enumerator$$Enumerator.prototype._enumerate = function() {
      var enumerator = this;

      var length  = enumerator.length;
      var promise = enumerator.promise;
      var input   = enumerator._input;

      for (var i = 0; promise._state === lib$es6$promise$$internal$$PENDING && i < length; i++) {
        enumerator._eachEntry(input[i], i);
      }
    };

    lib$es6$promise$enumerator$$Enumerator.prototype._eachEntry = function(entry, i) {
      var enumerator = this;
      var c = enumerator._instanceConstructor;

      if (lib$es6$promise$utils$$isMaybeThenable(entry)) {
        if (entry.constructor === c && entry._state !== lib$es6$promise$$internal$$PENDING) {
          entry._onerror = null;
          enumerator._settledAt(entry._state, i, entry._result);
        } else {
          enumerator._willSettleAt(c.resolve(entry), i);
        }
      } else {
        enumerator._remaining--;
        enumerator._result[i] = entry;
      }
    };

    lib$es6$promise$enumerator$$Enumerator.prototype._settledAt = function(state, i, value) {
      var enumerator = this;
      var promise = enumerator.promise;

      if (promise._state === lib$es6$promise$$internal$$PENDING) {
        enumerator._remaining--;

        if (state === lib$es6$promise$$internal$$REJECTED) {
          lib$es6$promise$$internal$$reject(promise, value);
        } else {
          enumerator._result[i] = value;
        }
      }

      if (enumerator._remaining === 0) {
        lib$es6$promise$$internal$$fulfill(promise, enumerator._result);
      }
    };

    lib$es6$promise$enumerator$$Enumerator.prototype._willSettleAt = function(promise, i) {
      var enumerator = this;

      lib$es6$promise$$internal$$subscribe(promise, undefined, function(value) {
        enumerator._settledAt(lib$es6$promise$$internal$$FULFILLED, i, value);
      }, function(reason) {
        enumerator._settledAt(lib$es6$promise$$internal$$REJECTED, i, reason);
      });
    };
    function lib$es6$promise$promise$all$$all(entries) {
      return new lib$es6$promise$enumerator$$default(this, entries).promise;
    }
    var lib$es6$promise$promise$all$$default = lib$es6$promise$promise$all$$all;
    function lib$es6$promise$promise$race$$race(entries) {
      /*jshint validthis:true */
      var Constructor = this;

      var promise = new Constructor(lib$es6$promise$$internal$$noop);

      if (!lib$es6$promise$utils$$isArray(entries)) {
        lib$es6$promise$$internal$$reject(promise, new TypeError('You must pass an array to race.'));
        return promise;
      }

      var length = entries.length;

      function onFulfillment(value) {
        lib$es6$promise$$internal$$resolve(promise, value);
      }

      function onRejection(reason) {
        lib$es6$promise$$internal$$reject(promise, reason);
      }

      for (var i = 0; promise._state === lib$es6$promise$$internal$$PENDING && i < length; i++) {
        lib$es6$promise$$internal$$subscribe(Constructor.resolve(entries[i]), undefined, onFulfillment, onRejection);
      }

      return promise;
    }
    var lib$es6$promise$promise$race$$default = lib$es6$promise$promise$race$$race;
    function lib$es6$promise$promise$resolve$$resolve(object) {
      /*jshint validthis:true */
      var Constructor = this;

      if (object && typeof object === 'object' && object.constructor === Constructor) {
        return object;
      }

      var promise = new Constructor(lib$es6$promise$$internal$$noop);
      lib$es6$promise$$internal$$resolve(promise, object);
      return promise;
    }
    var lib$es6$promise$promise$resolve$$default = lib$es6$promise$promise$resolve$$resolve;
    function lib$es6$promise$promise$reject$$reject(reason) {
      /*jshint validthis:true */
      var Constructor = this;
      var promise = new Constructor(lib$es6$promise$$internal$$noop);
      lib$es6$promise$$internal$$reject(promise, reason);
      return promise;
    }
    var lib$es6$promise$promise$reject$$default = lib$es6$promise$promise$reject$$reject;

    var lib$es6$promise$promise$$counter = 0;

    function lib$es6$promise$promise$$needsResolver() {
      throw new TypeError('You must pass a resolver function as the first argument to the promise constructor');
    }

    function lib$es6$promise$promise$$needsNew() {
      throw new TypeError("Failed to construct 'Promise': Please use the 'new' operator, this object constructor cannot be called as a function.");
    }

    var lib$es6$promise$promise$$default = lib$es6$promise$promise$$Promise;
    /**
      Promise objects represent the eventual result of an asynchronous operation. The
      primary way of interacting with a promise is through its `then` method, which
      registers callbacks to receive either a promise's eventual value or the reason
      why the promise cannot be fulfilled.

      Terminology
      -----------

      - `promise` is an object or function with a `then` method whose behavior conforms to this specification.
      - `thenable` is an object or function that defines a `then` method.
      - `value` is any legal JavaScript value (including undefined, a thenable, or a promise).
      - `exception` is a value that is thrown using the throw statement.
      - `reason` is a value that indicates why a promise was rejected.
      - `settled` the final resting state of a promise, fulfilled or rejected.

      A promise can be in one of three states: pending, fulfilled, or rejected.

      Promises that are fulfilled have a fulfillment value and are in the fulfilled
      state.  Promises that are rejected have a rejection reason and are in the
      rejected state.  A fulfillment value is never a thenable.

      Promises can also be said to *resolve* a value.  If this value is also a
      promise, then the original promise's settled state will match the value's
      settled state.  So a promise that *resolves* a promise that rejects will
      itself reject, and a promise that *resolves* a promise that fulfills will
      itself fulfill.


      Basic Usage:
      ------------

      ```js
      var promise = new Promise(function(resolve, reject) {
        // on success
        resolve(value);

        // on failure
        reject(reason);
      });

      promise.then(function(value) {
        // on fulfillment
      }, function(reason) {
        // on rejection
      });
      ```

      Advanced Usage:
      ---------------

      Promises shine when abstracting away asynchronous interactions such as
      `XMLHttpRequest`s.

      ```js
      function getJSON(url) {
        return new Promise(function(resolve, reject){
          var xhr = new XMLHttpRequest();

          xhr.open('GET', url);
          xhr.onreadystatechange = handler;
          xhr.responseType = 'json';
          xhr.setRequestHeader('Accept', 'application/json');
          xhr.send();

          function handler() {
            if (this.readyState === this.DONE) {
              if (this.status === 200) {
                resolve(this.response);
              } else {
                reject(new Error('getJSON: `' + url + '` failed with status: [' + this.status + ']'));
              }
            }
          };
        });
      }

      getJSON('/posts.json').then(function(json) {
        // on fulfillment
      }, function(reason) {
        // on rejection
      });
      ```

      Unlike callbacks, promises are great composable primitives.

      ```js
      Promise.all([
        getJSON('/posts'),
        getJSON('/comments')
      ]).then(function(values){
        values[0] // => postsJSON
        values[1] // => commentsJSON

        return values;
      });
      ```

      @class Promise
      @param {function} resolver
      Useful for tooling.
      @constructor
    */
    function lib$es6$promise$promise$$Promise(resolver) {
      this._id = lib$es6$promise$promise$$counter++;
      this._state = undefined;
      this._result = undefined;
      this._subscribers = [];

      if (lib$es6$promise$$internal$$noop !== resolver) {
        if (!lib$es6$promise$utils$$isFunction(resolver)) {
          lib$es6$promise$promise$$needsResolver();
        }

        if (!(this instanceof lib$es6$promise$promise$$Promise)) {
          lib$es6$promise$promise$$needsNew();
        }

        lib$es6$promise$$internal$$initializePromise(this, resolver);
      }
    }

    lib$es6$promise$promise$$Promise.all = lib$es6$promise$promise$all$$default;
    lib$es6$promise$promise$$Promise.race = lib$es6$promise$promise$race$$default;
    lib$es6$promise$promise$$Promise.resolve = lib$es6$promise$promise$resolve$$default;
    lib$es6$promise$promise$$Promise.reject = lib$es6$promise$promise$reject$$default;
    lib$es6$promise$promise$$Promise._setScheduler = lib$es6$promise$asap$$setScheduler;
    lib$es6$promise$promise$$Promise._setAsap = lib$es6$promise$asap$$setAsap;
    lib$es6$promise$promise$$Promise._asap = lib$es6$promise$asap$$asap;

    lib$es6$promise$promise$$Promise.prototype = {
      constructor: lib$es6$promise$promise$$Promise,

    /**
      The primary way of interacting with a promise is through its `then` method,
      which registers callbacks to receive either a promise's eventual value or the
      reason why the promise cannot be fulfilled.

      ```js
      findUser().then(function(user){
        // user is available
      }, function(reason){
        // user is unavailable, and you are given the reason why
      });
      ```

      Chaining
      --------

      The return value of `then` is itself a promise.  This second, 'downstream'
      promise is resolved with the return value of the first promise's fulfillment
      or rejection handler, or rejected if the handler throws an exception.

      ```js
      findUser().then(function (user) {
        return user.name;
      }, function (reason) {
        return 'default name';
      }).then(function (userName) {
        // If `findUser` fulfilled, `userName` will be the user's name, otherwise it
        // will be `'default name'`
      });

      findUser().then(function (user) {
        throw new Error('Found user, but still unhappy');
      }, function (reason) {
        throw new Error('`findUser` rejected and we're unhappy');
      }).then(function (value) {
        // never reached
      }, function (reason) {
        // if `findUser` fulfilled, `reason` will be 'Found user, but still unhappy'.
        // If `findUser` rejected, `reason` will be '`findUser` rejected and we're unhappy'.
      });
      ```
      If the downstream promise does not specify a rejection handler, rejection reasons will be propagated further downstream.

      ```js
      findUser().then(function (user) {
        throw new PedagogicalException('Upstream error');
      }).then(function (value) {
        // never reached
      }).then(function (value) {
        // never reached
      }, function (reason) {
        // The `PedgagocialException` is propagated all the way down to here
      });
      ```

      Assimilation
      ------------

      Sometimes the value you want to propagate to a downstream promise can only be
      retrieved asynchronously. This can be achieved by returning a promise in the
      fulfillment or rejection handler. The downstream promise will then be pending
      until the returned promise is settled. This is called *assimilation*.

      ```js
      findUser().then(function (user) {
        return findCommentsByAuthor(user);
      }).then(function (comments) {
        // The user's comments are now available
      });
      ```

      If the assimliated promise rejects, then the downstream promise will also reject.

      ```js
      findUser().then(function (user) {
        return findCommentsByAuthor(user);
      }).then(function (comments) {
        // If `findCommentsByAuthor` fulfills, we'll have the value here
      }, function (reason) {
        // If `findCommentsByAuthor` rejects, we'll have the reason here
      });
      ```

      Simple Example
      --------------

      Synchronous Example

      ```javascript
      var result;

      try {
        result = findResult();
        // success
      } catch(reason) {
        // failure
      }
      ```

      Errback Example

      ```js
      findResult(function(result, err){
        if (err) {
          // failure
        } else {
          // success
        }
      });
      ```

      Promise Example;

      ```javascript
      findResult().then(function(result){
        // success
      }, function(reason){
        // failure
      });
      ```

      Advanced Example
      --------------

      Synchronous Example

      ```javascript
      var author, books;

      try {
        author = findAuthor();
        books  = findBooksByAuthor(author);
        // success
      } catch(reason) {
        // failure
      }
      ```

      Errback Example

      ```js

      function foundBooks(books) {

      }

      function failure(reason) {

      }

      findAuthor(function(author, err){
        if (err) {
          failure(err);
          // failure
        } else {
          try {
            findBoooksByAuthor(author, function(books, err) {
              if (err) {
                failure(err);
              } else {
                try {
                  foundBooks(books);
                } catch(reason) {
                  failure(reason);
                }
              }
            });
          } catch(error) {
            failure(err);
          }
          // success
        }
      });
      ```

      Promise Example;

      ```javascript
      findAuthor().
        then(findBooksByAuthor).
        then(function(books){
          // found books
      }).catch(function(reason){
        // something went wrong
      });
      ```

      @method then
      @param {Function} onFulfilled
      @param {Function} onRejected
      Useful for tooling.
      @return {Promise}
    */
      then: function(onFulfillment, onRejection) {
        var parent = this;
        var state = parent._state;

        if (state === lib$es6$promise$$internal$$FULFILLED && !onFulfillment || state === lib$es6$promise$$internal$$REJECTED && !onRejection) {
          return this;
        }

        var child = new this.constructor(lib$es6$promise$$internal$$noop);
        var result = parent._result;

        if (state) {
          var callback = arguments[state - 1];
          lib$es6$promise$asap$$asap(function(){
            lib$es6$promise$$internal$$invokeCallback(state, child, callback, result);
          });
        } else {
          lib$es6$promise$$internal$$subscribe(parent, child, onFulfillment, onRejection);
        }

        return child;
      },

    /**
      `catch` is simply sugar for `then(undefined, onRejection)` which makes it the same
      as the catch block of a try/catch statement.

      ```js
      function findAuthor(){
        throw new Error('couldn't find that author');
      }

      // synchronous
      try {
        findAuthor();
      } catch(reason) {
        // something went wrong
      }

      // async with promises
      findAuthor().catch(function(reason){
        // something went wrong
      });
      ```

      @method catch
      @param {Function} onRejection
      Useful for tooling.
      @return {Promise}
    */
      'catch': function(onRejection) {
        return this.then(null, onRejection);
      }
    };
    function lib$es6$promise$polyfill$$polyfill() {
      var local;

      if (typeof global !== 'undefined') {
          local = global;
      } else if (typeof self !== 'undefined') {
          local = self;
      } else {
          try {
              local = Function('return this')();
          } catch (e) {
              throw new Error('polyfill failed because global object is unavailable in this environment');
          }
      }

      var P = local.Promise;

      if (P && Object.prototype.toString.call(P.resolve()) === '[object Promise]' && !P.cast) {
        return;
      }

      local.Promise = lib$es6$promise$promise$$default;
    }
    var lib$es6$promise$polyfill$$default = lib$es6$promise$polyfill$$polyfill;

    var lib$es6$promise$umd$$ES6Promise = {
      'Promise': lib$es6$promise$promise$$default,
      'polyfill': lib$es6$promise$polyfill$$default
    };

    /* global define:true module:true window: true */
    if (typeof define === 'function' && define['amd']) {
      define(function() { return lib$es6$promise$umd$$ES6Promise; });
    } else if (typeof module !== 'undefined' && module['exports']) {
      module['exports'] = lib$es6$promise$umd$$ES6Promise;
    } else if (typeof this !== 'undefined') {
      this['ES6Promise'] = lib$es6$promise$umd$$ES6Promise;
    }

    lib$es6$promise$polyfill$$default();
}).call(this);


}).call(this,require('_process'),typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})

},{"_process":14}],13:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Javascript ZLib
// By Thomas Down 2010-2011
//
// Based very heavily on portions of jzlib (by ymnk@jcraft.com), who in
// turn credits Jean-loup Gailly and Mark Adler for the original zlib code.
//
// inflate.js: ZLib inflate code
//

//
// Shared constants
//

var MAX_WBITS=15; // 32K LZ77 window
var DEF_WBITS=MAX_WBITS;
var MAX_MEM_LEVEL=9;
var MANY=1440;
var BMAX = 15;

// preset dictionary flag in zlib header
var PRESET_DICT=0x20;

var Z_NO_FLUSH=0;
var Z_PARTIAL_FLUSH=1;
var Z_SYNC_FLUSH=2;
var Z_FULL_FLUSH=3;
var Z_FINISH=4;

var Z_DEFLATED=8;

var Z_OK=0;
var Z_STREAM_END=1;
var Z_NEED_DICT=2;
var Z_ERRNO=-1;
var Z_STREAM_ERROR=-2;
var Z_DATA_ERROR=-3;
var Z_MEM_ERROR=-4;
var Z_BUF_ERROR=-5;
var Z_VERSION_ERROR=-6;

var METHOD=0;   // waiting for method byte
var FLAG=1;     // waiting for flag byte
var DICT4=2;    // four dictionary check bytes to go
var DICT3=3;    // three dictionary check bytes to go
var DICT2=4;    // two dictionary check bytes to go
var DICT1=5;    // one dictionary check byte to go
var DICT0=6;    // waiting for inflateSetDictionary
var BLOCKS=7;   // decompressing blocks
var CHECK4=8;   // four check bytes to go
var CHECK3=9;   // three check bytes to go
var CHECK2=10;  // two check bytes to go
var CHECK1=11;  // one check byte to go
var DONE=12;    // finished check, done
var BAD=13;     // got an error--stay here

var inflate_mask = [0x00000000, 0x00000001, 0x00000003, 0x00000007, 0x0000000f, 0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff, 0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff, 0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff];

var IB_TYPE=0;  // get type bits (3, including end bit)
var IB_LENS=1;  // get lengths for stored
var IB_STORED=2;// processing stored block
var IB_TABLE=3; // get table lengths
var IB_BTREE=4; // get bit lengths tree for a dynamic block
var IB_DTREE=5; // get length, distance trees for a dynamic block
var IB_CODES=6; // processing fixed or dynamic block
var IB_DRY=7;   // output remaining window bytes
var IB_DONE=8;  // finished last block, done
var IB_BAD=9;   // ot a data error--stuck here

var fixed_bl = 9;
var fixed_bd = 5;

var fixed_tl = [
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,192,
    80,7,10, 0,8,96, 0,8,32, 0,9,160,
    0,8,0, 0,8,128, 0,8,64, 0,9,224,
    80,7,6, 0,8,88, 0,8,24, 0,9,144,
    83,7,59, 0,8,120, 0,8,56, 0,9,208,
    81,7,17, 0,8,104, 0,8,40, 0,9,176,
    0,8,8, 0,8,136, 0,8,72, 0,9,240,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,200,
    81,7,13, 0,8,100, 0,8,36, 0,9,168,
    0,8,4, 0,8,132, 0,8,68, 0,9,232,
    80,7,8, 0,8,92, 0,8,28, 0,9,152,
    84,7,83, 0,8,124, 0,8,60, 0,9,216,
    82,7,23, 0,8,108, 0,8,44, 0,9,184,
    0,8,12, 0,8,140, 0,8,76, 0,9,248,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,196,
    81,7,11, 0,8,98, 0,8,34, 0,9,164,
    0,8,2, 0,8,130, 0,8,66, 0,9,228,
    80,7,7, 0,8,90, 0,8,26, 0,9,148,
    84,7,67, 0,8,122, 0,8,58, 0,9,212,
    82,7,19, 0,8,106, 0,8,42, 0,9,180,
    0,8,10, 0,8,138, 0,8,74, 0,9,244,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,204,
    81,7,15, 0,8,102, 0,8,38, 0,9,172,
    0,8,6, 0,8,134, 0,8,70, 0,9,236,
    80,7,9, 0,8,94, 0,8,30, 0,9,156,
    84,7,99, 0,8,126, 0,8,62, 0,9,220,
    82,7,27, 0,8,110, 0,8,46, 0,9,188,
    0,8,14, 0,8,142, 0,8,78, 0,9,252,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,194,
    80,7,10, 0,8,97, 0,8,33, 0,9,162,
    0,8,1, 0,8,129, 0,8,65, 0,9,226,
    80,7,6, 0,8,89, 0,8,25, 0,9,146,
    83,7,59, 0,8,121, 0,8,57, 0,9,210,
    81,7,17, 0,8,105, 0,8,41, 0,9,178,
    0,8,9, 0,8,137, 0,8,73, 0,9,242,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,202,
    81,7,13, 0,8,101, 0,8,37, 0,9,170,
    0,8,5, 0,8,133, 0,8,69, 0,9,234,
    80,7,8, 0,8,93, 0,8,29, 0,9,154,
    84,7,83, 0,8,125, 0,8,61, 0,9,218,
    82,7,23, 0,8,109, 0,8,45, 0,9,186,
    0,8,13, 0,8,141, 0,8,77, 0,9,250,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,198,
    81,7,11, 0,8,99, 0,8,35, 0,9,166,
    0,8,3, 0,8,131, 0,8,67, 0,9,230,
    80,7,7, 0,8,91, 0,8,27, 0,9,150,
    84,7,67, 0,8,123, 0,8,59, 0,9,214,
    82,7,19, 0,8,107, 0,8,43, 0,9,182,
    0,8,11, 0,8,139, 0,8,75, 0,9,246,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,206,
    81,7,15, 0,8,103, 0,8,39, 0,9,174,
    0,8,7, 0,8,135, 0,8,71, 0,9,238,
    80,7,9, 0,8,95, 0,8,31, 0,9,158,
    84,7,99, 0,8,127, 0,8,63, 0,9,222,
    82,7,27, 0,8,111, 0,8,47, 0,9,190,
    0,8,15, 0,8,143, 0,8,79, 0,9,254,
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,193,

    80,7,10, 0,8,96, 0,8,32, 0,9,161,
    0,8,0, 0,8,128, 0,8,64, 0,9,225,
    80,7,6, 0,8,88, 0,8,24, 0,9,145,
    83,7,59, 0,8,120, 0,8,56, 0,9,209,
    81,7,17, 0,8,104, 0,8,40, 0,9,177,
    0,8,8, 0,8,136, 0,8,72, 0,9,241,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,201,
    81,7,13, 0,8,100, 0,8,36, 0,9,169,
    0,8,4, 0,8,132, 0,8,68, 0,9,233,
    80,7,8, 0,8,92, 0,8,28, 0,9,153,
    84,7,83, 0,8,124, 0,8,60, 0,9,217,
    82,7,23, 0,8,108, 0,8,44, 0,9,185,
    0,8,12, 0,8,140, 0,8,76, 0,9,249,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,197,
    81,7,11, 0,8,98, 0,8,34, 0,9,165,
    0,8,2, 0,8,130, 0,8,66, 0,9,229,
    80,7,7, 0,8,90, 0,8,26, 0,9,149,
    84,7,67, 0,8,122, 0,8,58, 0,9,213,
    82,7,19, 0,8,106, 0,8,42, 0,9,181,
    0,8,10, 0,8,138, 0,8,74, 0,9,245,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,205,
    81,7,15, 0,8,102, 0,8,38, 0,9,173,
    0,8,6, 0,8,134, 0,8,70, 0,9,237,
    80,7,9, 0,8,94, 0,8,30, 0,9,157,
    84,7,99, 0,8,126, 0,8,62, 0,9,221,
    82,7,27, 0,8,110, 0,8,46, 0,9,189,
    0,8,14, 0,8,142, 0,8,78, 0,9,253,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,195,
    80,7,10, 0,8,97, 0,8,33, 0,9,163,
    0,8,1, 0,8,129, 0,8,65, 0,9,227,
    80,7,6, 0,8,89, 0,8,25, 0,9,147,
    83,7,59, 0,8,121, 0,8,57, 0,9,211,
    81,7,17, 0,8,105, 0,8,41, 0,9,179,
    0,8,9, 0,8,137, 0,8,73, 0,9,243,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,203,
    81,7,13, 0,8,101, 0,8,37, 0,9,171,
    0,8,5, 0,8,133, 0,8,69, 0,9,235,
    80,7,8, 0,8,93, 0,8,29, 0,9,155,
    84,7,83, 0,8,125, 0,8,61, 0,9,219,
    82,7,23, 0,8,109, 0,8,45, 0,9,187,
    0,8,13, 0,8,141, 0,8,77, 0,9,251,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,199,
    81,7,11, 0,8,99, 0,8,35, 0,9,167,
    0,8,3, 0,8,131, 0,8,67, 0,9,231,
    80,7,7, 0,8,91, 0,8,27, 0,9,151,
    84,7,67, 0,8,123, 0,8,59, 0,9,215,
    82,7,19, 0,8,107, 0,8,43, 0,9,183,
    0,8,11, 0,8,139, 0,8,75, 0,9,247,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,207,
    81,7,15, 0,8,103, 0,8,39, 0,9,175,
    0,8,7, 0,8,135, 0,8,71, 0,9,239,
    80,7,9, 0,8,95, 0,8,31, 0,9,159,
    84,7,99, 0,8,127, 0,8,63, 0,9,223,
    82,7,27, 0,8,111, 0,8,47, 0,9,191,
    0,8,15, 0,8,143, 0,8,79, 0,9,255
];
var fixed_td = [
    80,5,1, 87,5,257, 83,5,17, 91,5,4097,
    81,5,5, 89,5,1025, 85,5,65, 93,5,16385,
    80,5,3, 88,5,513, 84,5,33, 92,5,8193,
    82,5,9, 90,5,2049, 86,5,129, 192,5,24577,
    80,5,2, 87,5,385, 83,5,25, 91,5,6145,
    81,5,7, 89,5,1537, 85,5,97, 93,5,24577,
    80,5,4, 88,5,769, 84,5,49, 92,5,12289,
    82,5,13, 90,5,3073, 86,5,193, 192,5,24577
];

  // Tables for deflate from PKZIP's appnote.txt.
  var cplens = [ // Copy lengths for literal codes 257..285
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0
  ];

  // see note #13 above about 258
  var cplext = [ // Extra bits for literal codes 257..285
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 112, 112  // 112==invalid
  ];

 var cpdist = [ // Copy offsets for distance codes 0..29
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
  ];

  var cpdext = [ // Extra bits for distance codes
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13];

//
// ZStream.java
//

function ZStream() {
}


ZStream.prototype.inflateInit = function(w, nowrap) {
    if (!w) {
	w = DEF_WBITS;
    }
    if (nowrap) {
	nowrap = false;
    }
    this.istate = new Inflate();
    return this.istate.inflateInit(this, nowrap?-w:w);
}

ZStream.prototype.inflate = function(f) {
    if(this.istate==null) return Z_STREAM_ERROR;
    return this.istate.inflate(this, f);
}

ZStream.prototype.inflateEnd = function(){
    if(this.istate==null) return Z_STREAM_ERROR;
    var ret=istate.inflateEnd(this);
    this.istate = null;
    return ret;
}
ZStream.prototype.inflateSync = function(){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSync(this);
}
ZStream.prototype.inflateSetDictionary = function(dictionary, dictLength){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSetDictionary(this, dictionary, dictLength);
}

/*

  public int deflateInit(int level){
    return deflateInit(level, MAX_WBITS);
  }
  public int deflateInit(int level, boolean nowrap){
    return deflateInit(level, MAX_WBITS, nowrap);
  }
  public int deflateInit(int level, int bits){
    return deflateInit(level, bits, false);
  }
  public int deflateInit(int level, int bits, boolean nowrap){
    dstate=new Deflate();
    return dstate.deflateInit(this, level, nowrap?-bits:bits);
  }
  public int deflate(int flush){
    if(dstate==null){
      return Z_STREAM_ERROR;
    }
    return dstate.deflate(this, flush);
  }
  public int deflateEnd(){
    if(dstate==null) return Z_STREAM_ERROR;
    int ret=dstate.deflateEnd();
    dstate=null;
    return ret;
  }
  public int deflateParams(int level, int strategy){
    if(dstate==null) return Z_STREAM_ERROR;
    return dstate.deflateParams(this, level, strategy);
  }
  public int deflateSetDictionary (byte[] dictionary, int dictLength){
    if(dstate == null)
      return Z_STREAM_ERROR;
    return dstate.deflateSetDictionary(this, dictionary, dictLength);
  }

*/

/*
  // Flush as much pending output as possible. All deflate() output goes
  // through this function so some applications may wish to modify it
  // to avoid allocating a large strm->next_out buffer and copying into it.
  // (See also read_buf()).
  void flush_pending(){
    int len=dstate.pending;

    if(len>avail_out) len=avail_out;
    if(len==0) return;

    if(dstate.pending_buf.length<=dstate.pending_out ||
       next_out.length<=next_out_index ||
       dstate.pending_buf.length<(dstate.pending_out+len) ||
       next_out.length<(next_out_index+len)){
      System.out.println(dstate.pending_buf.length+", "+dstate.pending_out+
			 ", "+next_out.length+", "+next_out_index+", "+len);
      System.out.println("avail_out="+avail_out);
    }

    System.arraycopy(dstate.pending_buf, dstate.pending_out,
		     next_out, next_out_index, len);

    next_out_index+=len;
    dstate.pending_out+=len;
    total_out+=len;
    avail_out-=len;
    dstate.pending-=len;
    if(dstate.pending==0){
      dstate.pending_out=0;
    }
  }

  // Read a new buffer from the current input stream, update the adler32
  // and total number of bytes read.  All deflate() input goes through
  // this function so some applications may wish to modify it to avoid
  // allocating a large strm->next_in buffer and copying from it.
  // (See also flush_pending()).
  int read_buf(byte[] buf, int start, int size) {
    int len=avail_in;

    if(len>size) len=size;
    if(len==0) return 0;

    avail_in-=len;

    if(dstate.noheader==0) {
      adler=_adler.adler32(adler, next_in, next_in_index, len);
    }
    System.arraycopy(next_in, next_in_index, buf, start, len);
    next_in_index  += len;
    total_in += len;
    return len;
  }

  public void free(){
    next_in=null;
    next_out=null;
    msg=null;
    _adler=null;
  }
}
*/


//
// Inflate.java
//

function Inflate() {
    this.was = [0];
}

Inflate.prototype.inflateReset = function(z) {
    if(z == null || z.istate == null) return Z_STREAM_ERROR;
    
    z.total_in = z.total_out = 0;
    z.msg = null;
    z.istate.mode = z.istate.nowrap!=0 ? BLOCKS : METHOD;
    z.istate.blocks.reset(z, null);
    return Z_OK;
}

Inflate.prototype.inflateEnd = function(z){
    if(this.blocks != null)
      this.blocks.free(z);
    this.blocks=null;
    return Z_OK;
}

Inflate.prototype.inflateInit = function(z, w){
    z.msg = null;
    this.blocks = null;

    // handle undocumented nowrap option (no zlib header or check)
    nowrap = 0;
    if(w < 0){
      w = - w;
      nowrap = 1;
    }

    // set window size
    if(w<8 ||w>15){
      this.inflateEnd(z);
      return Z_STREAM_ERROR;
    }
    this.wbits=w;

    z.istate.blocks=new InfBlocks(z, 
				  z.istate.nowrap!=0 ? null : this,
				  1<<w);

    // reset state
    this.inflateReset(z);
    return Z_OK;
  }

Inflate.prototype.inflate = function(z, f){
    var r, b;

    if(z == null || z.istate == null || z.next_in == null)
      return Z_STREAM_ERROR;
    f = f == Z_FINISH ? Z_BUF_ERROR : Z_OK;
    r = Z_BUF_ERROR;
    while (true){
      switch (z.istate.mode){
      case METHOD:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        if(((z.istate.method = z.next_in[z.next_in_index++])&0xf)!=Z_DEFLATED){
          z.istate.mode = BAD;
          z.msg="unknown compression method";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        if((z.istate.method>>4)+8>z.istate.wbits){
          z.istate.mode = BAD;
          z.msg="invalid window size";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        z.istate.mode=FLAG;
      case FLAG:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        b = (z.next_in[z.next_in_index++])&0xff;

        if((((z.istate.method << 8)+b) % 31)!=0){
          z.istate.mode = BAD;
          z.msg = "incorrect header check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        if((b&PRESET_DICT)==0){
          z.istate.mode = BLOCKS;
          break;
        }
        z.istate.mode = DICT4;
      case DICT4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=DICT3;
      case DICT3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode=DICT2;
      case DICT2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode=DICT1;
      case DICT1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need += (z.next_in[z.next_in_index++]&0xff);
        z.adler = z.istate.need;
        z.istate.mode = DICT0;
        return Z_NEED_DICT;
      case DICT0:
        z.istate.mode = BAD;
        z.msg = "need dictionary";
        z.istate.marker = 0;       // can try inflateSync
        return Z_STREAM_ERROR;
      case BLOCKS:

        r = z.istate.blocks.proc(z, r);
        if(r == Z_DATA_ERROR){
          z.istate.mode = BAD;
          z.istate.marker = 0;     // can try inflateSync
          break;
        }
        if(r == Z_OK){
          r = f;
        }
        if(r != Z_STREAM_END){
          return r;
        }
        r = f;
        z.istate.blocks.reset(z, z.istate.was);
        if(z.istate.nowrap!=0){
          z.istate.mode=DONE;
          break;
        }
        z.istate.mode=CHECK4;
      case CHECK4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=CHECK3;
      case CHECK3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode = CHECK2;
      case CHECK2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode = CHECK1;
      case CHECK1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=(z.next_in[z.next_in_index++]&0xff);

        if(((z.istate.was[0])) != ((z.istate.need))){
          z.istate.mode = BAD;
          z.msg = "incorrect data check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        z.istate.mode = DONE;
      case DONE:
        return Z_STREAM_END;
      case BAD:
        return Z_DATA_ERROR;
      default:
        return Z_STREAM_ERROR;
      }
    }
  }


Inflate.prototype.inflateSetDictionary = function(z,  dictionary, dictLength) {
    var index=0;
    var length = dictLength;
    if(z==null || z.istate == null|| z.istate.mode != DICT0)
      return Z_STREAM_ERROR;

    if(z._adler.adler32(1, dictionary, 0, dictLength)!=z.adler){
      return Z_DATA_ERROR;
    }

    z.adler = z._adler.adler32(0, null, 0, 0);

    if(length >= (1<<z.istate.wbits)){
      length = (1<<z.istate.wbits)-1;
      index=dictLength - length;
    }
    z.istate.blocks.set_dictionary(dictionary, index, length);
    z.istate.mode = BLOCKS;
    return Z_OK;
  }

//  static private byte[] mark = {(byte)0, (byte)0, (byte)0xff, (byte)0xff};
var mark = [0, 0, 255, 255]

Inflate.prototype.inflateSync = function(z){
    var n;       // number of bytes to look at
    var p;       // pointer to bytes
    var m;       // number of marker bytes found in a row
    var r, w;   // temporaries to save total_in and total_out

    // set up
    if(z == null || z.istate == null)
      return Z_STREAM_ERROR;
    if(z.istate.mode != BAD){
      z.istate.mode = BAD;
      z.istate.marker = 0;
    }
    if((n=z.avail_in)==0)
      return Z_BUF_ERROR;
    p=z.next_in_index;
    m=z.istate.marker;

    // search
    while (n!=0 && m < 4){
      if(z.next_in[p] == mark[m]){
        m++;
      }
      else if(z.next_in[p]!=0){
        m = 0;
      }
      else{
        m = 4 - m;
      }
      p++; n--;
    }

    // restore
    z.total_in += p-z.next_in_index;
    z.next_in_index = p;
    z.avail_in = n;
    z.istate.marker = m;

    // return no joy or set up to restart on a new block
    if(m != 4){
      return Z_DATA_ERROR;
    }
    r=z.total_in;  w=z.total_out;
    this.inflateReset(z);
    z.total_in=r;  z.total_out = w;
    z.istate.mode = BLOCKS;
    return Z_OK;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. This function is used by one PPP
  // implementation to provide an additional safety check. PPP uses Z_SYNC_FLUSH
  // but removes the length bytes of the resulting empty stored block. When
  // decompressing, PPP checks that at the end of input packet, inflate is
  // waiting for these length bytes.
Inflate.prototype.inflateSyncPoint = function(z){
    if(z == null || z.istate == null || z.istate.blocks == null)
      return Z_STREAM_ERROR;
    return z.istate.blocks.sync_point();
}


//
// InfBlocks.java
//

var INFBLOCKS_BORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

function InfBlocks(z, checkfn, w) {
    this.hufts=new Int32Array(MANY*3);
    this.window=new Uint8Array(w);
    this.end=w;
    this.checkfn = checkfn;
    this.mode = IB_TYPE;
    this.reset(z, null);

    this.left = 0;            // if STORED, bytes left to copy 

    this.table = 0;           // table lengths (14 bits) 
    this.index = 0;           // index into blens (or border) 
    this.blens = null;         // bit lengths of codes 
    this.bb=new Int32Array(1); // bit length tree depth 
    this.tb=new Int32Array(1); // bit length decoding tree 

    this.codes = new InfCodes();

    this.last = 0;            // true if this block is the last block 

  // mode independent information 
    this.bitk = 0;            // bits in bit buffer 
    this.bitb = 0;            // bit buffer 
    this.read = 0;            // window read pointer 
    this.write = 0;           // window write pointer 
    this.check = 0;          // check on output 

    this.inftree=new InfTree();
}




InfBlocks.prototype.reset = function(z, c){
    if(c) c[0]=this.check;
    if(this.mode==IB_CODES){
      this.codes.free(z);
    }
    this.mode=IB_TYPE;
    this.bitk=0;
    this.bitb=0;
    this.read=this.write=0;

    if(this.checkfn)
      z.adler=this.check=z._adler.adler32(0, null, 0, 0);
  }

 InfBlocks.prototype.proc = function(z, r){
    var t;              // temporary storage
    var b;              // bit buffer
    var k;              // bits in bit buffer
    var p;              // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer

    // copy input/output information to locals (UPDATE macro restores)
    {p=z.next_in_index;n=z.avail_in;b=this.bitb;k=this.bitk;}
    {q=this.write;m=(q<this.read ? this.read-q-1 : this.end-q);}

    // process input based on current state
    while(true){
      switch (this.mode){
      case IB_TYPE:

	while(k<(3)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}
	t = (b & 7);
	this.last = t & 1;

	switch (t >>> 1){
        case 0:                         // stored 
          {b>>>=(3);k-=(3);}
          t = k & 7;                    // go to byte boundary

          {b>>>=(t);k-=(t);}
          this.mode = IB_LENS;                  // get length of stored block
          break;
        case 1:                         // fixed
          {
              var bl=new Int32Array(1);
	      var bd=new Int32Array(1);
              var tl=[];
	      var td=[];

	      inflate_trees_fixed(bl, bd, tl, td, z);
              this.codes.init(bl[0], bd[0], tl[0], 0, td[0], 0, z);
          }

          {b>>>=(3);k-=(3);}

          this.mode = IB_CODES;
          break;
        case 2:                         // dynamic

          {b>>>=(3);k-=(3);}

          this.mode = IB_TABLE;
          break;
        case 3:                         // illegal

          {b>>>=(3);k-=(3);}
          this.mode = BAD;
          z.msg = "invalid block type";
          r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	break;
      case IB_LENS:
	while(k<(32)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	if ((((~b) >>> 16) & 0xffff) != (b & 0xffff)){
	  this.mode = BAD;
	  z.msg = "invalid stored block lengths";
	  r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	this.left = (b & 0xffff);
	b = k = 0;                       // dump bits
	this.mode = this.left!=0 ? IB_STORED : (this.last!=0 ? IB_DRY : IB_TYPE);
	break;
      case IB_STORED:
	if (n == 0){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	if(m==0){
	  if(q==end&&read!=0){
	    q=0; m=(q<this.read ? this.read-q-1 : this.end-q);
	  }
	  if(m==0){
	    this.write=q; 
	    r=this.inflate_flush(z,r);
	    q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	    if(q==this.end && this.read != 0){
	      q=0; m = (q < this.read ? this.read-q-1 : this.end-q);
	    }
	    if(m==0){
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	t = this.left;
	if(t>n) t = n;
	if(t>m) t = m;
	arrayCopy(z.next_in, p, this.window, q, t);
	p += t;  n -= t;
	q += t;  m -= t;
	if ((this.left -= t) != 0)
	  break;
	this.mode = (this.last != 0 ? IB_DRY : IB_TYPE);
	break;
      case IB_TABLE:

	while(k<(14)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.table = t = (b & 0x3fff);
	if ((t & 0x1f) > 29 || ((t >> 5) & 0x1f) > 29)
	  {
	    this.mode = IB_BAD;
	    z.msg = "too many length or distance symbols";
	    r = Z_DATA_ERROR;

	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  }
	t = 258 + (t & 0x1f) + ((t >> 5) & 0x1f);
	if(this.blens==null || this.blens.length<t){
	    this.blens=new Int32Array(t);
	}
	else{
	  for(var i=0; i<t; i++){
              this.blens[i]=0;
          }
	}

	{b>>>=(14);k-=(14);}

	this.index = 0;
	mode = IB_BTREE;
      case IB_BTREE:
	while (this.index < 4 + (this.table >>> 10)){
	  while(k<(3)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

	  this.blens[INFBLOCKS_BORDER[this.index++]] = b&7;

	  {b>>>=(3);k-=(3);}
	}

	while(this.index < 19){
	  this.blens[INFBLOCKS_BORDER[this.index++]] = 0;
	}

	this.bb[0] = 7;
	t = this.inftree.inflate_trees_bits(this.blens, this.bb, this.tb, this.hufts, z);
	if (t != Z_OK){
	  r = t;
	  if (r == Z_DATA_ERROR){
	    this.blens=null;
	    this.mode = IB_BAD;
	  }

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	this.index = 0;
	this.mode = IB_DTREE;
      case IB_DTREE:
	while (true){
	  t = this.table;
	  if(!(this.index < 258 + (t & 0x1f) + ((t >> 5) & 0x1f))){
	    break;
	  }

	  var h; //int[]
	  var i, j, c;

	  t = this.bb[0];

	  while(k<(t)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

//	  if (this.tb[0]==-1){
//            dlog("null...");
//	  }

	  t=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+1];
	  c=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+2];

	  if (c < 16){
	    b>>>=(t);k-=(t);
	    this.blens[this.index++] = c;
	  }
	  else { // c == 16..18
	    i = c == 18 ? 7 : c - 14;
	    j = c == 18 ? 11 : 3;

	    while(k<(t+i)){
	      if(n!=0){
		r=Z_OK;
	      }
	      else{
		this.bitb=b; this.bitk=k; 
		z.avail_in=n;
		z.total_in+=p-z.next_in_index;z.next_in_index=p;
		this.write=q;
		return this.inflate_flush(z,r);
	      };
	      n--;
	      b|=(z.next_in[p++]&0xff)<<k;
	      k+=8;
	    }

	    b>>>=(t);k-=(t);

	    j += (b & inflate_mask[i]);

	    b>>>=(i);k-=(i);

	    i = this.index;
	    t = this.table;
	    if (i + j > 258 + (t & 0x1f) + ((t >> 5) & 0x1f) ||
		(c == 16 && i < 1)){
	      this.blens=null;
	      this.mode = IB_BAD;
	      z.msg = "invalid bit length repeat";
	      r = Z_DATA_ERROR;

	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }

	    c = c == 16 ? this.blens[i-1] : 0;
	    do{
	      this.blens[i++] = c;
	    }
	    while (--j!=0);
	    this.index = i;
	  }
	}

	this.tb[0]=-1;
	{
	    var bl=new Int32Array(1);
	    var bd=new Int32Array(1);
	    var tl=new Int32Array(1);
	    var td=new Int32Array(1);
	    bl[0] = 9;         // must be <= 9 for lookahead assumptions
	    bd[0] = 6;         // must be <= 9 for lookahead assumptions

	    t = this.table;
	    t = this.inftree.inflate_trees_dynamic(257 + (t & 0x1f), 
					      1 + ((t >> 5) & 0x1f),
					      this.blens, bl, bd, tl, td, this.hufts, z);

	    if (t != Z_OK){
	        if (t == Z_DATA_ERROR){
	            this.blens=null;
	            this.mode = BAD;
	        }
	        r = t;

	        this.bitb=b; this.bitk=k; 
	        z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	        this.write=q;
	        return this.inflate_flush(z,r);
	    }
	    this.codes.init(bl[0], bd[0], this.hufts, tl[0], this.hufts, td[0], z);
	}
	this.mode = IB_CODES;
      case IB_CODES:
	this.bitb=b; this.bitk=k;
	z.avail_in=n; z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;

	if ((r = this.codes.proc(this, z, r)) != Z_STREAM_END){
	  return this.inflate_flush(z, r);
	}
	r = Z_OK;
	this.codes.free(z);

	p=z.next_in_index; n=z.avail_in;b=this.bitb;k=this.bitk;
	q=this.write;m = (q < this.read ? this.read-q-1 : this.end-q);

	if (this.last==0){
	  this.mode = IB_TYPE;
	  break;
	}
	this.mode = IB_DRY;
      case IB_DRY:
	this.write=q; 
	r = this.inflate_flush(z, r); 
	q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	if (this.read != this.write){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z, r);
	}
	mode = DONE;
      case IB_DONE:
	r = Z_STREAM_END;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      case IB_BAD:
	r = Z_DATA_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);

      default:
	r = Z_STREAM_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      }
    }
  }

InfBlocks.prototype.free = function(z){
    this.reset(z, null);
    this.window=null;
    this.hufts=null;
}

InfBlocks.prototype.set_dictionary = function(d, start, n){
    arrayCopy(d, start, window, 0, n);
    this.read = this.write = n;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. 
InfBlocks.prototype.sync_point = function(){
    return this.mode == IB_LENS;
}

  // copy as much as possible from the sliding window to the output area
InfBlocks.prototype.inflate_flush = function(z, r){
    var n;
    var p;
    var q;

    // local copies of source and destination pointers
    p = z.next_out_index;
    q = this.read;

    // compute number of bytes to copy as far as end of window
    n = ((q <= this.write ? this.write : this.end) - q);
    if (n > z.avail_out) n = z.avail_out;
    if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

    // update counters
    z.avail_out -= n;
    z.total_out += n;

    // update check information
    if(this.checkfn != null)
      z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

    // copy as far as end of window
    arrayCopy(this.window, q, z.next_out, p, n);
    p += n;
    q += n;

    // see if more to copy at beginning of window
    if (q == this.end){
      // wrap pointers
      q = 0;
      if (this.write == this.end)
        this.write = 0;

      // compute bytes to copy
      n = this.write - q;
      if (n > z.avail_out) n = z.avail_out;
      if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

      // update counters
      z.avail_out -= n;
      z.total_out += n;

      // update check information
      if(this.checkfn != null)
	z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

      // copy
      arrayCopy(this.window, q, z.next_out, p, n);
      p += n;
      q += n;
    }

    // update pointers
    z.next_out_index = p;
    this.read = q;

    // done
    return r;
  }

//
// InfCodes.java
//

var IC_START=0;  // x: set up for LEN
var IC_LEN=1;    // i: get length/literal/eob next
var IC_LENEXT=2; // i: getting length extra (have base)
var IC_DIST=3;   // i: get distance next
var IC_DISTEXT=4;// i: getting distance extra
var IC_COPY=5;   // o: copying bytes in window, waiting for space
var IC_LIT=6;    // o: got literal, waiting for output space
var IC_WASH=7;   // o: got eob, possibly still output waiting
var IC_END=8;    // x: got eob and all data flushed
var IC_BADCODE=9;// x: got error

function InfCodes() {
}

InfCodes.prototype.init = function(bl, bd, tl, tl_index, td, td_index, z) {
    this.mode=IC_START;
    this.lbits=bl;
    this.dbits=bd;
    this.ltree=tl;
    this.ltree_index=tl_index;
    this.dtree = td;
    this.dtree_index=td_index;
    this.tree=null;
}

InfCodes.prototype.proc = function(s, z, r){ 
    var j;              // temporary storage
    var t;              // temporary pointer (int[])
    var tindex;         // temporary pointer
    var e;              // extra bits or operation
    var b=0;            // bit buffer
    var k=0;            // bits in bit buffer
    var p=0;            // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer
    var f;              // pointer to copy strings from

    // copy input/output information to locals (UPDATE macro restores)
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // process input and output based on current state
    while (true){
      switch (this.mode){
	// waiting for "i:"=input, "o:"=output, "x:"=nothing
      case IC_START:         // x: set up for LEN
	if (m >= 258 && n >= 10){

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  r = this.inflate_fast(this.lbits, this.dbits, 
			   this.ltree, this.ltree_index, 
			   this.dtree, this.dtree_index,
			   s, z);

	  p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
	  q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	  if (r != Z_OK){
	    this.mode = r == Z_STREAM_END ? IC_WASH : IC_BADCODE;
	    break;
	  }
	}
	this.need = this.lbits;
	this.tree = this.ltree;
	this.tree_index=this.ltree_index;

	this.mode = IC_LEN;
      case IC_LEN:           // i: get length/literal/eob next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b&inflate_mask[j]))*3;

	b>>>=(this.tree[tindex+1]);
	k-=(this.tree[tindex+1]);

	e=this.tree[tindex];

	if(e == 0){               // literal
	  this.lit = this.tree[tindex+2];
	  this.mode = IC_LIT;
	  break;
	}
	if((e & 16)!=0 ){          // length
	  this.get = e & 15;
	  this.len = this.tree[tindex+2];
	  this.mode = IC_LENEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	if ((e & 32)!=0){               // end of block
	  this.mode = IC_WASH;
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid literal/length code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_LENEXT:        // i: getting length extra (have base)
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.len += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.need = this.dbits;
	this.tree = this.dtree;
	this.tree_index = this.dtree_index;
	this.mode = IC_DIST;
      case IC_DIST:          // i: get distance next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b & inflate_mask[j]))*3;

	b>>=this.tree[tindex+1];
	k-=this.tree[tindex+1];

	e = (this.tree[tindex]);
	if((e & 16)!=0){               // distance
	  this.get = e & 15;
	  this.dist = this.tree[tindex+2];
	  this.mode = IC_DISTEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid distance code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_DISTEXT:       // i: getting distance extra
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.dist += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.mode = IC_COPY;
      case IC_COPY:          // o: copying bytes in window, waiting for space
        f = q - this.dist;
        while(f < 0){     // modulo window size-"while" instead
          f += s.end;     // of "if" handles invalid distances
	}
	while (this.len!=0){

	  if(m==0){
	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.write=q; r=s.inflate_flush(z,r);
	      q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	      if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}

	      if(m==0){
		s.bitb=b;s.bitk=k;
		z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
		s.write=q;
		return s.inflate_flush(z,r);
	      }  
	    }
	  }

	  s.window[q++]=s.window[f++]; m--;

	  if (f == s.end)
            f = 0;
	  this.len--;
	}
	this.mode = IC_START;
	break;
      case IC_LIT:           // o: got literal, waiting for output space
	if(m==0){
	  if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	  if(m==0){
	    s.write=q; r=s.inflate_flush(z,r);
	    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;
	      return s.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	s.window[q++]=this.lit; m--;

	this.mode = IC_START;
	break;
      case IC_WASH:           // o: got eob, possibly more output
	if (k > 7){        // return unused byte, if any
	  k -= 8;
	  n++;
	  p--;             // can always return one
	}

	s.write=q; r=s.inflate_flush(z,r);
	q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	if (s.read != s.write){
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  return s.inflate_flush(z,r);
	}
	this.mode = IC_END;
      case IC_END:
	r = Z_STREAM_END;
	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_BADCODE:       // x: got error

	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      default:
	r = Z_STREAM_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);
      }
    }
  }

InfCodes.prototype.free = function(z){
    //  ZFREE(z, c);
}

  // Called with number of bytes left to write in window at least 258
  // (the maximum string length) and number of input bytes available
  // at least ten.  The ten bytes are six bytes for the longest length/
  // distance pair plus four bytes for overloading the bit buffer.

InfCodes.prototype.inflate_fast = function(bl, bd, tl, tl_index, td, td_index, s, z) {
    var t;                // temporary pointer
    var   tp;             // temporary pointer (int[])
    var tp_index;         // temporary pointer
    var e;                // extra bits or operation
    var b;                // bit buffer
    var k;                // bits in bit buffer
    var p;                // input data pointer
    var n;                // bytes available there
    var q;                // output window write pointer
    var m;                // bytes to end of window or read pointer
    var ml;               // mask for literal/length tree
    var md;               // mask for distance tree
    var c;                // bytes to copy
    var d;                // distance back to copy from
    var r;                // copy source pointer

    var tp_index_t_3;     // (tp_index+t)*3

    // load input, output, bit values
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // initialize masks
    ml = inflate_mask[bl];
    md = inflate_mask[bd];

    // do until not enough input or output space for fast loop
    do {                          // assume called with m >= 258 && n >= 10
      // get literal/length code
      while(k<(20)){              // max bits for literal/length code
	n--;
	b|=(z.next_in[p++]&0xff)<<k;k+=8;
      }

      t= b&ml;
      tp=tl; 
      tp_index=tl_index;
      tp_index_t_3=(tp_index+t)*3;
      if ((e = tp[tp_index_t_3]) == 0){
	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	s.window[q++] = tp[tp_index_t_3+2];
	m--;
	continue;
      }
      do {

	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	if((e&16)!=0){
	  e &= 15;
	  c = tp[tp_index_t_3+2] + (b & inflate_mask[e]);

	  b>>=e; k-=e;

	  // decode distance base of block to copy
	  while(k<(15)){           // max bits for distance code
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;k+=8;
	  }

	  t= b&md;
	  tp=td;
	  tp_index=td_index;
          tp_index_t_3=(tp_index+t)*3;
	  e = tp[tp_index_t_3];

	  do {

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    if((e&16)!=0){
	      // get extra bits to add to distance base
	      e &= 15;
	      while(k<(e)){         // get extra bits (up to 13)
		n--;
		b|=(z.next_in[p++]&0xff)<<k;k+=8;
	      }

	      d = tp[tp_index_t_3+2] + (b&inflate_mask[e]);

	      b>>=(e); k-=(e);

	      // do the copy
	      m -= c;
	      if (q >= d){                // offset before dest
		//  just copy
		r=q-d;
		if(q-r>0 && 2>(q-r)){           
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
		else{
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
	      }
	      else{                  // else offset after destination
                r=q-d;
                do{
                  r+=s.end;          // force pointer in window
                }while(r<0);         // covers invalid distances
		e=s.end-r;
		if(c>e){             // if source crosses,
		  c-=e;              // wrapped copy
		  if(q-r>0 && e>(q-r)){           
		    do{s.window[q++] = s.window[r++];}
		    while(--e!=0);
		  }
		  else{
		    arrayCopy(s.window, r, s.window, q, e);
		    q+=e; r+=e; e=0;
		  }
		  r = 0;                  // copy rest from start of window
		}

	      }

	      // copy all or what's left
              do{s.window[q++] = s.window[r++];}
		while(--c!=0);
	      break;
	    }
	    else if((e&64)==0){
	      t+=tp[tp_index_t_3+2];
	      t+=(b&inflate_mask[e]);
	      tp_index_t_3=(tp_index+t)*3;
	      e=tp[tp_index_t_3];
	    }
	    else{
	      z.msg = "invalid distance code";

	      c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;

	      return Z_DATA_ERROR;
	    }
	  }
	  while(true);
	  break;
	}

	if((e&64)==0){
	  t+=tp[tp_index_t_3+2];
	  t+=(b&inflate_mask[e]);
	  tp_index_t_3=(tp_index+t)*3;
	  if((e=tp[tp_index_t_3])==0){

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    s.window[q++]=tp[tp_index_t_3+2];
	    m--;
	    break;
	  }
	}
	else if((e&32)!=0){

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;
 
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_STREAM_END;
	}
	else{
	  z.msg="invalid literal/length code";

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_DATA_ERROR;
	}
      } 
      while(true);
    } 
    while(m>=258 && n>= 10);

    // not enough input or output--restore pointers and return
    c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

    s.bitb=b;s.bitk=k;
    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
    s.write=q;

    return Z_OK;
}

//
// InfTree.java
//

function InfTree() {
}

InfTree.prototype.huft_build = function(b, bindex, n, s, d, e, t, m, hp, hn, v) {

    // Given a list of code lengths and a maximum table size, make a set of
    // tables to decode that set of codes.  Return Z_OK on success, Z_BUF_ERROR
    // if the given code set is incomplete (the tables are still built in this
    // case), Z_DATA_ERROR if the input is invalid (an over-subscribed set of
    // lengths), or Z_MEM_ERROR if not enough memory.

    var a;                       // counter for codes of length k
    var f;                       // i repeats in table every f entries
    var g;                       // maximum code length
    var h;                       // table level
    var i;                       // counter, current code
    var j;                       // counter
    var k;                       // number of bits in current code
    var l;                       // bits per table (returned in m)
    var mask;                    // (1 << w) - 1, to avoid cc -O bug on HP
    var p;                       // pointer into c[], b[], or v[]
    var q;                       // points to current table
    var w;                       // bits before this table == (l * h)
    var xp;                      // pointer into x
    var y;                       // number of dummy codes added
    var z;                       // number of entries in current table

    // Generate counts for each bit length

    p = 0; i = n;
    do {
      this.c[b[bindex+p]]++; p++; i--;   // assume all entries <= BMAX
    }while(i!=0);

    if(this.c[0] == n){                // null input--all zero length codes
      t[0] = -1;
      m[0] = 0;
      return Z_OK;
    }

    // Find minimum and maximum length, bound *m by those
    l = m[0];
    for (j = 1; j <= BMAX; j++)
      if(this.c[j]!=0) break;
    k = j;                        // minimum code length
    if(l < j){
      l = j;
    }
    for (i = BMAX; i!=0; i--){
      if(this.c[i]!=0) break;
    }
    g = i;                        // maximum code length
    if(l > i){
      l = i;
    }
    m[0] = l;

    // Adjust last length count to fill out codes, if needed
    for (y = 1 << j; j < i; j++, y <<= 1){
      if ((y -= this.c[j]) < 0){
        return Z_DATA_ERROR;
      }
    }
    if ((y -= this.c[i]) < 0){
      return Z_DATA_ERROR;
    }
    this.c[i] += y;

    // Generate starting offsets into the value table for each length
    this.x[1] = j = 0;
    p = 1;  xp = 2;
    while (--i!=0) {                 // note that i == g from above
      this.x[xp] = (j += this.c[p]);
      xp++;
      p++;
    }

    // Make a table of values in order of bit lengths
    i = 0; p = 0;
    do {
      if ((j = b[bindex+p]) != 0){
        this.v[this.x[j]++] = i;
      }
      p++;
    }
    while (++i < n);
    n = this.x[g];                     // set n to length of v

    // Generate the Huffman codes and for each, make the table entries
    this.x[0] = i = 0;                 // first Huffman code is zero
    p = 0;                        // grab values in bit order
    h = -1;                       // no tables yet--level -1
    w = -l;                       // bits decoded == (l * h)
    this.u[0] = 0;                     // just to keep compilers happy
    q = 0;                        // ditto
    z = 0;                        // ditto

    // go through the bit lengths (k already is bits in shortest code)
    for (; k <= g; k++){
      a = this.c[k];
      while (a--!=0){
	// here i is the Huffman code of length k bits for value *p
	// make tables up to required level
        while (k > w + l){
          h++;
          w += l;                 // previous table always l bits
	  // compute minimum size table less than or equal to l bits
          z = g - w;
          z = (z > l) ? l : z;        // table size upper limit
          if((f=1<<(j=k-w))>a+1){     // try a k-w bit table
                                      // too few codes for k-w bit table
            f -= a + 1;               // deduct codes from patterns left
            xp = k;
            if(j < z){
              while (++j < z){        // try smaller tables up to z bits
                if((f <<= 1) <= this.c[++xp])
                  break;              // enough codes to use up j bits
                f -= this.c[xp];           // else deduct codes from patterns
              }
	    }
          }
          z = 1 << j;                 // table entries for j-bit table

	  // allocate new table
          if (this.hn[0] + z > MANY){       // (note: doesn't matter for fixed)
            return Z_DATA_ERROR;       // overflow of MANY
          }
          this.u[h] = q = /*hp+*/ this.hn[0];   // DEBUG
          this.hn[0] += z;
 
	  // connect to last table, if there is one
	  if(h!=0){
            this.x[h]=i;           // save pattern for backing up
            this.r[0]=j;     // bits in this table
            this.r[1]=l;     // bits to dump before this table
            j=i>>>(w - l);
            this.r[2] = (q - this.u[h-1] - j);               // offset to this table
            arrayCopy(this.r, 0, hp, (this.u[h-1]+j)*3, 3); // connect to last table
          }
          else{
            t[0] = q;               // first table is returned result
	  }
        }

	// set up table entry in r
        this.r[1] = (k - w);
        if (p >= n){
          this.r[0] = 128 + 64;      // out of values--invalid code
	}
        else if (v[p] < s){
          this.r[0] = (this.v[p] < 256 ? 0 : 32 + 64);  // 256 is end-of-block
          this.r[2] = this.v[p++];          // simple code is just the value
        }
        else{
          this.r[0]=(e[this.v[p]-s]+16+64); // non-simple--look up in lists
          this.r[2]=d[this.v[p++] - s];
        }

        // fill code-like entries with r
        f=1<<(k-w);
        for (j=i>>>w;j<z;j+=f){
          arrayCopy(this.r, 0, hp, (q+j)*3, 3);
	}

	// backwards increment the k-bit code i
        for (j = 1 << (k - 1); (i & j)!=0; j >>>= 1){
          i ^= j;
	}
        i ^= j;

	// backup over finished tables
        mask = (1 << w) - 1;      // needed on HP, cc -O bug
        while ((i & mask) != this.x[h]){
          h--;                    // don't need to update q
          w -= l;
          mask = (1 << w) - 1;
        }
      }
    }
    // Return Z_BUF_ERROR if we were given an incomplete table
    return y != 0 && g != 1 ? Z_BUF_ERROR : Z_OK;
}

InfTree.prototype.inflate_trees_bits = function(c, bb, tb, hp, z) {
    var result;
    this.initWorkArea(19);
    this.hn[0]=0;
    result = this.huft_build(c, 0, 19, 19, null, null, tb, bb, hp, this.hn, this.v);

    if(result == Z_DATA_ERROR){
      z.msg = "oversubscribed dynamic bit lengths tree";
    }
    else if(result == Z_BUF_ERROR || bb[0] == 0){
      z.msg = "incomplete dynamic bit lengths tree";
      result = Z_DATA_ERROR;
    }
    return result;
}

InfTree.prototype.inflate_trees_dynamic = function(nl, nd, c, bl, bd, tl, td, hp, z) {
    var result;

    // build literal/length tree
    this.initWorkArea(288);
    this.hn[0]=0;
    result = this.huft_build(c, 0, nl, 257, cplens, cplext, tl, bl, hp, this.hn, this.v);
    if (result != Z_OK || bl[0] == 0){
      if(result == Z_DATA_ERROR){
        z.msg = "oversubscribed literal/length tree";
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "incomplete literal/length tree";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    // build distance tree
    this.initWorkArea(288);
    result = this.huft_build(c, nl, nd, 0, cpdist, cpdext, td, bd, hp, this.hn, this.v);

    if (result != Z_OK || (bd[0] == 0 && nl > 257)){
      if (result == Z_DATA_ERROR){
        z.msg = "oversubscribed distance tree";
      }
      else if (result == Z_BUF_ERROR) {
        z.msg = "incomplete distance tree";
        result = Z_DATA_ERROR;
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "empty distance tree with lengths";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    return Z_OK;
}
/*
  static int inflate_trees_fixed(int[] bl,  //literal desired/actual bit depth
                                 int[] bd,  //distance desired/actual bit depth
                                 int[][] tl,//literal/length tree result
                                 int[][] td,//distance tree result 
                                 ZStream z  //for memory allocation
				 ){

*/

function inflate_trees_fixed(bl, bd, tl, td, z) {
    bl[0]=fixed_bl;
    bd[0]=fixed_bd;
    tl[0]=fixed_tl;
    td[0]=fixed_td;
    return Z_OK;
}

InfTree.prototype.initWorkArea = function(vsize){
    if(this.hn==null){
        this.hn=new Int32Array(1);
        this.v=new Int32Array(vsize);
        this.c=new Int32Array(BMAX+1);
        this.r=new Int32Array(3);
        this.u=new Int32Array(BMAX);
        this.x=new Int32Array(BMAX+1);
    }
    if(this.v.length<vsize){ 
        this.v=new Int32Array(vsize); 
    }
    for(var i=0; i<vsize; i++){this.v[i]=0;}
    for(var i=0; i<BMAX+1; i++){this.c[i]=0;}
    for(var i=0; i<3; i++){this.r[i]=0;}
//  for(int i=0; i<BMAX; i++){u[i]=0;}
    arrayCopy(this.c, 0, this.u, 0, BMAX);
//  for(int i=0; i<BMAX+1; i++){x[i]=0;}
    arrayCopy(this.c, 0, this.x, 0, BMAX+1);
}

var testArray = new Uint8Array(1);
var hasSubarray = (typeof testArray.subarray === 'function');
var hasSlice = false; /* (typeof testArray.slice === 'function'); */ // Chrome slice performance is so dire that we're currently not using it...

function arrayCopy(src, srcOffset, dest, destOffset, count) {
    if (count == 0) {
        return;
    } 
    if (!src) {
        throw "Undef src";
    } else if (!dest) {
        throw "Undef dest";
    }

    if (srcOffset == 0 && count == src.length) {
        arrayCopy_fast(src, dest, destOffset);
    } else if (hasSubarray) {
        arrayCopy_fast(src.subarray(srcOffset, srcOffset + count), dest, destOffset); 
    } else if (src.BYTES_PER_ELEMENT == 1 && count > 100) {
        arrayCopy_fast(new Uint8Array(src.buffer, src.byteOffset + srcOffset, count), dest, destOffset);
    } else { 
        arrayCopy_slow(src, srcOffset, dest, destOffset, count);
    }

}

function arrayCopy_slow(src, srcOffset, dest, destOffset, count) {

    // dlog('_slow call: srcOffset=' + srcOffset + '; destOffset=' + destOffset + '; count=' + count);

     for (var i = 0; i < count; ++i) {
        dest[destOffset + i] = src[srcOffset + i];
    }
}

function arrayCopy_fast(src, dest, destOffset) {
    dest.set(src, destOffset);
}


  // largest prime smaller than 65536
var ADLER_BASE=65521; 
  // NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1
var ADLER_NMAX=5552;

function adler32(adler, /* byte[] */ buf,  index, len){
    if(buf == null){ return 1; }

    var s1=adler&0xffff;
    var s2=(adler>>16)&0xffff;
    var k;

    while(len > 0) {
      k=len<ADLER_NMAX?len:ADLER_NMAX;
      len-=k;
      while(k>=16){
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        k-=16;
      }
      if(k!=0){
        do{
          s1+=buf[index++]&0xff; s2+=s1;
        }
        while(--k!=0);
      }
      s1%=ADLER_BASE;
      s2%=ADLER_BASE;
    }
    return (s2<<16)|s1;
}



function jszlib_inflate_buffer(buffer, start, length, afterUncOffset) {
    if (!start) {
        buffer = new Uint8Array(buffer);
    } else if (!length) {
        buffer = new Uint8Array(buffer, start, buffer.byteLength - start);
    } else {
        buffer = new Uint8Array(buffer, start, length);
    }

    var z = new ZStream();
    z.inflateInit(DEF_WBITS, true);
    z.next_in = buffer;
    z.next_in_index = 0;
    z.avail_in = buffer.length;

    var oBlockList = [];
    var totalSize = 0;
    while (true) {
        var obuf = new Uint8Array(32000);
        z.next_out = obuf;
        z.next_out_index = 0;
        z.avail_out = obuf.length;
        var status = z.inflate(Z_NO_FLUSH);
        if (status != Z_OK && status != Z_STREAM_END && status != Z_BUF_ERROR) {
            throw z.msg;
        }
        if (z.avail_out != 0) {
            var newob = new Uint8Array(obuf.length - z.avail_out);
            arrayCopy(obuf, 0, newob, 0, (obuf.length - z.avail_out));
            obuf = newob;
        }
        oBlockList.push(obuf);
        totalSize += obuf.length;
        if (status == Z_STREAM_END || status == Z_BUF_ERROR) {
            break;
        }
    }

    if (afterUncOffset) {
        afterUncOffset[0] = (start || 0) + z.next_in_index;
    }

    if (oBlockList.length == 1) {
        return oBlockList[0].buffer;
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = oBlockList[i];
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

if (typeof(module) !== 'undefined') {
  module.exports = {
    inflateBuffer: jszlib_inflate_buffer,
    arrayCopy: arrayCopy
  };
}

},{}],14:[function(require,module,exports){
// shim for using process in browser
var process = module.exports = {};

// cached from whatever global is present so that test runners that stub it
// don't break things.  But we need to wrap it in a try catch in case it is
// wrapped in strict mode code which doesn't define any globals.  It's inside a
// function because try/catches deoptimize in certain engines.

var cachedSetTimeout;
var cachedClearTimeout;

function defaultSetTimout() {
    throw new Error('setTimeout has not been defined');
}
function defaultClearTimeout () {
    throw new Error('clearTimeout has not been defined');
}
(function () {
    try {
        if (typeof setTimeout === 'function') {
            cachedSetTimeout = setTimeout;
        } else {
            cachedSetTimeout = defaultSetTimout;
        }
    } catch (e) {
        cachedSetTimeout = defaultSetTimout;
    }
    try {
        if (typeof clearTimeout === 'function') {
            cachedClearTimeout = clearTimeout;
        } else {
            cachedClearTimeout = defaultClearTimeout;
        }
    } catch (e) {
        cachedClearTimeout = defaultClearTimeout;
    }
} ())
function runTimeout(fun) {
    if (cachedSetTimeout === setTimeout) {
        //normal enviroments in sane situations
        return setTimeout(fun, 0);
    }
    // if setTimeout wasn't available but was latter defined
    if ((cachedSetTimeout === defaultSetTimout || !cachedSetTimeout) && setTimeout) {
        cachedSetTimeout = setTimeout;
        return setTimeout(fun, 0);
    }
    try {
        // when when somebody has screwed with setTimeout but no I.E. maddness
        return cachedSetTimeout(fun, 0);
    } catch(e){
        try {
            // When we are in I.E. but the script has been evaled so I.E. doesn't trust the global object when called normally
            return cachedSetTimeout.call(null, fun, 0);
        } catch(e){
            // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error
            return cachedSetTimeout.call(this, fun, 0);
        }
    }


}
function runClearTimeout(marker) {
    if (cachedClearTimeout === clearTimeout) {
        //normal enviroments in sane situations
        return clearTimeout(marker);
    }
    // if clearTimeout wasn't available but was latter defined
    if ((cachedClearTimeout === defaultClearTimeout || !cachedClearTimeout) && clearTimeout) {
        cachedClearTimeout = clearTimeout;
        return clearTimeout(marker);
    }
    try {
        // when when somebody has screwed with setTimeout but no I.E. maddness
        return cachedClearTimeout(marker);
    } catch (e){
        try {
            // When we are in I.E. but the script has been evaled so I.E. doesn't  trust the global object when called normally
            return cachedClearTimeout.call(null, marker);
        } catch (e){
            // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error.
            // Some versions of I.E. have different rules for clearTimeout vs setTimeout
            return cachedClearTimeout.call(this, marker);
        }
    }



}
var queue = [];
var draining = false;
var currentQueue;
var queueIndex = -1;

function cleanUpNextTick() {
    if (!draining || !currentQueue) {
        return;
    }
    draining = false;
    if (currentQueue.length) {
        queue = currentQueue.concat(queue);
    } else {
        queueIndex = -1;
    }
    if (queue.length) {
        drainQueue();
    }
}

function drainQueue() {
    if (draining) {
        return;
    }
    var timeout = runTimeout(cleanUpNextTick);
    draining = true;

    var len = queue.length;
    while(len) {
        currentQueue = queue;
        queue = [];
        while (++queueIndex < len) {
            if (currentQueue) {
                currentQueue[queueIndex].run();
            }
        }
        queueIndex = -1;
        len = queue.length;
    }
    currentQueue = null;
    draining = false;
    runClearTimeout(timeout);
}

process.nextTick = function (fun) {
    var args = new Array(arguments.length - 1);
    if (arguments.length > 1) {
        for (var i = 1; i < arguments.length; i++) {
            args[i - 1] = arguments[i];
        }
    }
    queue.push(new Item(fun, args));
    if (queue.length === 1 && !draining) {
        runTimeout(drainQueue);
    }
};

// v8 likes predictible objects
function Item(fun, array) {
    this.fun = fun;
    this.array = array;
}
Item.prototype.run = function () {
    this.fun.apply(null, this.array);
};
process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];
process.version = ''; // empty string to avoid regexp issues
process.versions = {};

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;
process.prependListener = noop;
process.prependOnceListener = noop;

process.listeners = function (name) { return [] }

process.binding = function (name) {
    throw new Error('process.binding is not supported');
};

process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};
process.umask = function() { return 0; };

},{}]},{},[7])

//# sourceMappingURL=worker-all.js.map
