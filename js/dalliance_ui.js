
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
      hg19ToHg38: new Chainset('http://www.derkholm.net:8080/das/hg19ToHg38/', 'GRCh37', 'GRCh38',
                    {
                        speciesName: 'Human',
                        taxon: 9606,
                        auth: 'GRCh',
                        version: 37,
                        ucscName: 'hg19'
                    })
    },

    browserLinks: {
      Ensembl: 'http://www.ensembl.org/Homo_sapiens/Location/View?r=${chr}:${start}-${end}',
      UCSC: 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${chr}:${start}-${end}',
      Sequence: 'http://www.derkholm.net:8080/das/hg19ToHg38/sequence?segment=${chr}:${start},${end}'
    },

    hubs: ['http://ngs.sanger.ac.uk/production/ensembl/regulation/hub.txt', {url: 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt', genome: 'hg19', mapping: 'hg19ToHg38'}],

	sources: [
            {name: 'Genome',
            twoBitURI:            'http://www.biodalliance.org/datasets/hg38.2bit',
            tier_type: 'sequence'},

            {name: 'Repeats',
            disabled: true,
            desc: 'Repeat annotation from UCSC', 
            bwgURI: '//www.biodalliance.org/datasets/GRCh38/repeats.bb',
            stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats2.xml'},

            {name: 'GENCODEv21',
            desc: 'Gene structures from GENCODE 21',
            bwgURI: '//www.biodalliance.org/datasets/GRCh38/gencode.v21.annotation.bb',
            stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode2.xml',
            collapseSuperGroups: true,
            trixURI: '//www.biodalliance.org/datasets/GRCh38/gencode.v21.annotation.ix'},

            {name: 'Ensembl REST',
            desc: 'Ensembl Transcript Model of GRCh38',
            uri: 'http://rest.ensembl.org',
            collapseSuperGroups: true,
            tier_type: 'ensembl',
            species: 'human',
            type: ['gene','transcript','exon','cds']},

            {name: 'CHESS v2.1',
            desc: 'CHESS v2.1',
            collapseSuperGroups: true,
            bwgURI: '//nectarmutation.org/main/static/chess2.1_assembly.ens.sorted.bb',
            stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode2.xml'},

            {name: 'POPS Transcriptome v82',
            desc: 'Placenta Transcriptome',
            collapseSuperGroups: true,
            bwgURI: '//nectarmutation.org/main/static/POPS.GRCh38.82.cuffcompare.combined.bb',
            stylesheet_uri: 'http://www.biodalliance.org/stylesheets/gencode.xml',
            trixURI: '//nectarmutation.org/main/static/POPS.GRCh38.82.cuffcompare.combined.geneName.ix'},

            {name: 'POPS30 circRNA',
            desc: 'circular RNAs present in 30% of POPS cohort',
            bwgURI: '//nectarmutation.org/main/static/circRNA.POPS30.CIRI2.GRCh38.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'darkgreen', BGCOLOR: 'rgb(117,214,146)', LABEL:'yes'}}]},

            {name: 'POPS novel miRNA',
            desc: 'novel miRNAs present in at least one sample of POPS cohort',
            bwgURI: '//nectarmutation.org/main/static/novel.miRNA.GRCh38.merged.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'orange', BGCOLOR: 'yellow', LABEL:'yes'}}]},

            {name: 'POPS novel small-RNA',
            desc: 'novel small-RNA present in at least one sample of POPS cohort',
            bwgURI: '//nectarmutation.org/main/static/novel.smallRNA.GRCh38.merged.bb',
            style: [{type: 'default', style:{glyph:'ANCHORED_ARROW', FGCOLOR: 'navy', BGCOLOR: 'blue', LABEL:'yes'}}]},

            {name: 'POPS total-RNA Coverage',
            desc: 'Total RNA-Seq Coverage of POPS Cohort',
            bwgURI: '//nectarmutation.org/main/static/POPS.GRCh38.82.bw',
            style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(166,71,71)', HEIGHT: 60, id: 'style1'}}],
            noDownsample: true,
            },

	],
  });

