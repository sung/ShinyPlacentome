
  var b = new Browser({
    chr:                 '19',
    viewStart:           666365,
    viewEnd:             693399,
    cookieKey:           'human-grc_h38',

    coordSystem: {
      speciesName: 'Human',
      taxon: 9606,
      auth: 'GRCh',
      version: '38',
      ucscName: 'hg38'
    },

    uiPrefix: '//www.biodalliance.org/release-0.13/',
    fullScreen: true,
    reverseScrolling: true,

    chains: {
        hg19ToHg38: {
            // uri: '//www.biodalliance.org/datasets/hg19ToHg38.bb',
            uri: '//www.obgyn.cam.ac.uk/Biodalliance/hg19ToHg38.bb',
            type: 'bigbed',
            coords: {
                speciesName: 'Human',
                taxon: 9606,
                auth: 'GRCh',
                version: 37,
                ucscName: 'hg19'
            }
        }
    },

    browserLinks: {
      Ensembl: '//www.ensembl.org/Homo_sapiens/Location/View?r=${chr}:${start}-${end}',
      UCSC: '//genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${chr}:${start}-${end}',
      Sequence: '//www.derkholm.net:8080/das/hg19ToHg38/sequence?segment=${chr}:${start},${end}'
    },

    hubs: [
            {url: '//ftp.ensembl.org/pub/papers/regulation/hub.txt', 
            forceProtocol: 'https',
            genome: 'hg38'
            //mapping: 'hg19ToHg38'
            },

            {url: 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt',
            forceProtocol: 'http',
            genome: 'hg19',
            mapping: 'hg19ToHg38'}
    ],

	sources: [
            {name: 'Genome',
            twoBitURI: '//www.obgyn.cam.ac.uk/Biodalliance/GRCh38.2bit',
            tier_type: 'sequence'},

            {name: 'Repeats',
            disabled: true,
            desc: 'Repeat annotation from UCSC', 
            //bwgURI: '//www.biodalliance.org/datasets/GRCh38/repeats.bb',
            //stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats2.xml'},
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/repeats.bb',
            stylesheet_uri: '//www.obgyn.cam.ac.uk/Biodalliance/bb-repeats2.xml'},

            {name: 'GENCODEv21',
            disabled: true,
            desc: 'Gene structures from GENCODE 21',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/gencode.v21.annotation.bb',
            stylesheet_uri: '//www.obgyn.cam.ac.uk/Biodalliance/gencode2.xml',
            collapseSuperGroups: true,
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/gencode.v21.annotation.ix'},

            {name: 'Ensembl REST',
            disabled: true,
            desc: 'Ensembl Transcript Model of GRCh38',
            uri: '//rest.ensembl.org',
            collapseSuperGroups: true,
            tier_type: 'ensembl',
            species: 'human',
            type: ['gene','transcript','exon','cds']},
        
            {name: 'Ensembl v94',
            desc: 'Ensembl Transcript Model v94',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/Homo_sapiens.GRCh38.94.bb',
            stylesheet_uri: '//www.obgyn.cam.ac.uk/Biodalliance/gencode2.xml',
            collapseSuperGroups: true,
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/Homo_sapiens.GRCh38.94.geneName.ix'},

            {name: 'CHESS v2.1',
            desc: 'CHESS v2.1',
            collapseSuperGroups: true,
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/chess2.1_assembly.ens.sorted.bb',
            stylesheet_uri: '//www.obgyn.cam.ac.uk/Biodalliance/gencode2.xml',
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/chess2.1_assembly.ens.sorted.geneName.ix'},

            {name: 'POPS Transcriptome v82',
            desc: 'Placenta Transcriptome',
            collapseSuperGroups: true,
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/POPS.GRCh38.82.cuffcompare.combined.geneName.bb',
            stylesheet_uri: '//www.obgyn.cam.ac.uk/Biodalliance/gencode.xml',
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/POPS.GRCh38.82.cuffcompare.combined.geneName.ix'},

            {name: 'POPS30 circRNA',
            desc: 'circular RNAs present in 30% of POPS cohort',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/circRNA.POPS30.CIRI2.GRCh38.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'darkgreen', BGCOLOR: 'rgb(117,214,146)', LABEL:'yes',HEIGHT:7}}]},

            {name: 'POPS total-RNA Coverage',
            desc: 'Total RNA-Seq Coverage of POPS Cohort',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/POPS.GRCh38.82.bw',
            style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(166,71,71)', HEIGHT: 60, id: 'style1'}}],
            noDownsample: true,
            },

            {name: 'miRBase v21',
            desc: 'microRNA sequences and annotations from miRBase v21',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/miRbase.v21.GRCh38.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'rgb(196,86,39)', BGCOLOR: 'rgb(255,214,104)', LABEL:'yes',HEIGHT:8}}],
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/miRbase.v21.GRCh38.ix'},

            {name: 'miRBase v22',
            desc: 'microRNA sequences and annotations from miRBase v22',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/miRbase.22.GRCh38.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'rgb(204,98,53)', BGCOLOR: 'rgb(255,214,104)', LABEL:'yes',HEIGHT:8}}]},

            {name: 'piRBase v1.0',
            desc: 'piwi-interacting RNA sequences and annotations from piRBase v1',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/piRNA.v1-miR.overlap30.GRCh38.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'navy', BGCOLOR: 'rgb(13, 187, 221)', LABEL:'yes',HEIGHT:8}}],
            trixURI: '//www.obgyn.cam.ac.uk/Biodalliance/piRNA.v1-miR.overlap30.GRCh38.ID.ix'},
        
            {name: 'POPS novel miRNA',
            desc: 'novel miRNAs present in at least one sample of POPS cohort',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/novel.miRNA.GRCh38.merged.at.least.two.methods.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'orange', BGCOLOR: 'yellow', LABEL:'yes',HEIGHT:8}}]},

            {name: 'POPS novel small-RNA',
            desc: 'novel small-RNA present in at least one sample of POPS cohort',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/novel.smallRNA.GRCh38.merged.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'navy', BGCOLOR: 'blue', LABEL:'yes',HEIGHT:8}}]},

            {name: 'POPS small-RNA Coverage',
            desc: 'small RNA-Seq Coverage of POPS Cohort',
            bwgURI: '//www.obgyn.cam.ac.uk/Biodalliance/POPS.small-RNA.GRCh38.depth.more.than.0.bw',
            style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(75,71,166)', HEIGHT: 60, id: 'style1'}}],
            noDownsample: true}
	],
  });

