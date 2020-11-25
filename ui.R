# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019

library(shiny)
library(shinythemes)
library(DT)
library(markdown) # The ‘includeMarkdown’ function requires the ‘markdown’ package
library(d3heatmap)

#load("RData/DEG.RData") # client-side gene selection with selectize
#gene.names<-dt.deseq[!is.na(hgnc_symbol),.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol # client side list of genes

navbarPage(title="POPS Placenta Transcriptome",
    theme=shinytheme("sandstone"),
    # below works for 'dalliance.js'
#    withTags({
#        head(
#                script(
#                    src="http://www.biodalliance.org/release-0.13/dalliance-compiled.js"
#                )
#        )
#    }),
#    tag$head(tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js")), # does not work

    tabPanel("Home",
        includeMarkdown("home.md")
    ),

    navbarMenu("Browse DEG",
        tabPanel("By P-value",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Select differentially expressed genes by a statistical significance level (p-value)"),
                        # drop down - transcript type
                        selectInput("transcript_pval", 
                                label="Choose a type of transcript:",
                                choices = list(
                                            "protein coding"="total-RNA:protein_coding", 
                                            "lincRNA"="total-RNA:lincRNA",
                                            "circRNA"="circRNA",
                                            "miRNA"="miRNA",
                                            "piRNA"="piRNA",
                                            "small non-coding RNA (snc-RNA)"="sncRNA",
                                            "novel protein-coding RNA"="novel-RPT:protein_coding",
                                            "novel lincRNA"="novel-RPT:lincRNA",
                                            "novel miRNA"="novel_miRNA",
                                            "novel small-RNA"="novel_smallRNA"),
                                selected="total-RNA:protein_coding"
                        ),
                        # drop down - pvalue 
                        selectInput("pval", 
                                    label = "Adjusted P-value:", 
                                    choices=list(`<0.1`=0.1,`<0.05`=0.05,`<0.01`=0.01),
                                    selected=0.05),
                        # Download button
                        downloadButton("download_pval", "Download")
                    ), # end of sidebarPanel
                    # Show a plot of the generated distribution
                    mainPanel(
                        tabsetPanel(
                            tabPanel("PE",
                                     DT::dataTableOutput('deg_pval_PE')
                            ),
                            tabPanel("SGA",
                                     DT::dataTableOutput('deg_pval_SGA')
                            )
                        )
                    ) # end of mainPanel
                ) # end of sidebarLayout
            ) # end of fluidPage
        ), # end of tabPanel p-value

        tabPanel("By fold-change (bootstrap)",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Select differentially expressed genes within top 5% fold-change calculated by bootstrapping of samples 5000 times"),
                        # drop down 
                        selectInput("transcript_boot", 
                                label="Choose a type of transcript:",
                                choices = list(
                                            "protein coding"="total-RNA:protein_coding", 
                                            "lincRNA"="total-RNA:lincRNA",
                                            "circRNA"="circRNA",
                                            "miRNA"="miRNA",
                                            "piRNA"="piRNA",
                                            "small non-coding RNA (snc-RNA)"="sncRNA",
                                            "novel protein-coding RNA"="novel-RPT:protein_coding",
                                            "novel lincRNA"="novel-RPT:lincRNA",
                                            "novel miRNA"="novel_miRNA",
                                            "novel small-RNA"="novel_smallRNA"),
                                selected="total-RNA:protein_coding"
                        ),
                        # Download button
                        downloadButton("download_boot", "Download")
                    ), # end of sidebarPanel
                    # Show a plot of the generated distribution
                    mainPanel(
                        tabsetPanel(
                            tabPanel("PE",
                                     #DT::dataTableOutput('deg_boot_pe_all'),
                                     #hr(),
                                     DT::dataTableOutput('deg_boot_pe_oneThird')
                            ),
                            tabPanel("SGA",
                                     #DT::dataTableOutput('deg_boot_sga_all'),
                                     #hr(),
                                     DT::dataTableOutput('deg_boot_sga_oneThird')
                            )
                        )
                    ) # end of mainPanel
                ) # end of sidebarLayout
            ) # end of fluidPage
        ), # end of tabPanel by FC

        tabPanel("By gene names",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Search DEGs by gene names(s) HGNC"),
                    # drop down 
                    #    textInput("gene", "Gene Name (HGNC):", "FSTL3") # old-school
                    #    selectizeInput("deg_genes","Gene Name(s):", choices=gene.names, selected = "FSTL4", multiple=TRUE) # client side
                        selectizeInput("deg_genes","Gene Name(s):", choices=NULL, selected="FSTL3", multiple=TRUE) # server side
                    ),
                    # Show a plot of the generated distribution
                    mainPanel(
                        tabsetPanel(
                            tabPanel("By P-value",
                                DT::dataTableOutput('deg_gene_pval')
                            ),
                            tabPanel("By Bootstrap",
                                DT::dataTableOutput('deg_gene_boot')
                            )
                        )
                    ) # end of mainPanel
                )
            )
        ), # end of tabPanel - by gene name

        tabPanel("By ENSEMBL IDs",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Search DEGs by Ensembl ID(s)"),
                        # drop down 
                        selectizeInput("deg_ensgs","ENSG ID(s):", choices=NULL, multiple=TRUE) # server side
                    ),
                    # Show a plot of the generated distribution
                    mainPanel(
                        tabsetPanel(
                            tabPanel("By P-value",
                                DT::dataTableOutput('deg_ensg_pval')
                            ),
                            tabPanel("By Bootstrap",
                                DT::dataTableOutput('deg_ensg_boot')
                            )
                        )
                    ) # end of mainPanel
                )
            )
        ) # end of tabPanel by Ensembl
    ), # end of navbarMenu DEG
    
    tabPanel(title="Browse Transcript",
        fluidPage(
            sidebarLayout(
                sidebarPanel(
                    helpText("Browse the transcript abundance level"),
                    #helpText("Browse the abundance level of potentially novel re-constructed transcripts"),
                    # drop down 
                    selectInput("ab_transcript", 
                            label="Choose a type of transcript:",
                            choices = list(
                                        "protein coding"="protein_coding",  # total-RNA
                                        "lincRNA"="lincRNA",
                                        "pseudogene"="pseudogene",
                                        "circRNA"="circRNA",
                                        "miRNA"="miRNA",
                                        "piRNA"="piRNA",
                                        "small non-coding RNA (snc-RNA)"="sncRNA",
                                        "novel transcript isoforms"="novel_isoform",
                                        "novel miRNA"="novel_miRNA",
                                        "novel small-RNA"="novel_smallRNA"),
                            selected="protein_coding"),
                    # conditional drop down - min FPKM of placenta (except circRNA)
                    conditionalPanel(
                        condition="input.ab_transcript != 'circRNA'",
                        selectInput("fpkm", 
                                    label = "Minimum FPKM:", 
                                    choices=list(`>=0`=0,`>0.01`=0.01,`>0.1`=0.1,`>1`=1,`>5`=5,`>10`=10),
                                    selected=1)),
                    # conditional select button - polyA+ or not (circRNA only)
                    conditionalPanel(
                        condition="input.ab_transcript == 'circRNA'",
                        checkboxInput("in_polyA", label = "EXCLUDE circRNA found in polyA+ data?", value = TRUE)),
                    # conditional no. of exon - only for novel_isoform
                    conditionalPanel(
                        condition="input.ab_transcript == 'novel_isoform'",
                        checkboxInput("no_single_exon", label = "EXCLUDE single-exon transcript?", value = TRUE)),
                    # conditional sample frequency range (only for circRNA, novel_isoform)
                    conditionalPanel(
                        condition="input.ab_transcript == 'novel_isoform' | input.ab_transcript =='circRNA'",
                        sliderInput("evi.ratio", label = "Transcript present in following sample frequency range", min = 0, max = 1, value = c(0.3,1))),
                    downloadButton("download_pops_tr", "Download")
                    # a set of slide
                ),
                # Show a plot of the generated distribution
                mainPanel(
                    DT::dataTableOutput('pops_tr')
                ) # end of mainPanel
            ) # end of sidebarLayout
        ) # end of fludPage
    ), # end of tabPanel - Browse Transcript

    navbarMenu("Placenta vs GTEx",
        tabPanel(title="Placenta enriched",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Browse transcripts enriched in the placenta"),
                        # drop down 
                        selectInput("transcript_tau", 
                                label="Choose a type of transcript (Ensembl v90 or Gencode v27):",
                                choices = list(
                                            "protein coding"="protein_coding", 
                                            "lincRNA"="lincRNA",
                                            "processed pseudogene"="processed_pseudogene")),
                        # drop down - min FPKM of placenta 
                        selectInput("pt_tpm1", 
                                    label = "Minimum TPM of Placenta:", 
                                    choices=list(`>=0.1`=0.1,`>=1`=1,`>=5`=5, `>=10`=10),
                                    selected=1),
                        # a set of radio buttons - transcript type
                        # tau score range 
                        sliderInput("tau", label = a("Tau score",href="https://academic.oup.com/bib/article/18/2/205/2562739",target="_blank"), min = 0.5, max = 1, value = c(0.99,1)),
                        # drop down - min FPKM of placenta 
                        selectInput("pt_gtex_fc", 
                                    label = "Fold change of placenta compared with the average of 49 GTEx tissues:", 
                                    choices=list(`>10x`=10,`>100x`=100,`>1000x`=1000),
                                    selected=100),
                        downloadButton("download_tau", "Download")
                        # a set of radio buttons - transcript type
                    ),
                    # Show a plot of the generated distribution
                    mainPanel(
                        #verbatimTextOutput("options"),
                        #verbatimTextOutput("test4"),
                        verbatimTextOutput("heatmap_title"),
                        d3heatmapOutput("heatmap", width="100%", height="1100px"),
                        hr(),
                        DT::dataTableOutput('tau')
                    ) # end of mainPanel
                ) # end of sidebarLayout
            ) # end of fludPage
        ), # end of tabPanel placenta-specific

        tabPanel("genes of your interests",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        helpText("Tissue-wide comparision of expression level"),
                        # radio button
                        radioButtons("radio_gtex", label ="List all or search genes of your interests",
                            choices = list("All" = 1, "Gene(s) of your interests" = 2), 
                            selected = 1),
                        # drop down - min FPKM of placenta 
                        conditionalPanel(
                            condition="input.radio_gtex==1",
                            selectInput("pt_tpm2", 
                                        label = "Minimum TPM of Placenta:", 
                                        choices=list(`>=0`=0,`>=0.1`=0.1,`>=1`=1,`>=5`=5, `>=10`=10,`>=100`=100),
                                        selected=1)
                        ),
                        conditionalPanel(
                            condition="input.radio_gtex==2",
                            selectizeInput("gtex_genes","Gene Name(s):", choices=NULL, selected="FSTL3", multiple=TRUE) # server side
                        ),
                        downloadButton("download_gtex", "Download")
                    ),
                    # Show a plot of the generated distribution
                    mainPanel(
                        DT::dataTableOutput('gtex_tpm'),
                        conditionalPanel(
                            condition="input.radio_gtex==2",
                            plotOutput("gtex_tpm_barchart")
                        )
                    ) # end of mainPanel
                )
            )
        ) # end of tabPanel - by gene name
#        ), # end of tabPanel - by gene name
#
#        tabPanel("Not in placenta",
#            fluidPage(
#                sidebarLayout(
#                    sidebarPanel(
#                        helpText("Browse genes not enriched in the placneta compared other tissues"),
#                        # checkbox of tissues to exclude
#                        checkboxGroupInput("no_gtex_tissue", 
#                                      label = "EXCLUDE following tissues from GTEx", 
#                                      choices=gtex_tissues,
#                                      selected=c("Blood","Breast")
#                                      ),
#                        # drop down 
#                        selectInput("transcript_not_in_pt", 
#                                label="Choose a type of transcript (Ensembl v90 or Gencode v27):",
#                                choices = list("protein coding"="protein_coding", "lincRNA"="lincRNA"),
#                                selected="protein_coding"),
#                        # checkbox
#                        conditionalPanel(
#                            condition="input.transcript_not_in_pt== 'protein_coding'",
#                            checkboxInput("no_ribosomal", label = "EXCLUDE ribosomal protein?", value = FALSE)),
#                        # drop down - min baseMean
#                        selectInput("min_gtex_count", 
#                                    label = "Minimum read count of non-placental tissue:", 
#                                    choices=list(`>10`=10,`>20`=20,`>50`=50, `>100`=100),
#                                    selected=10),
#                        # drop down - min TPM of GTEx 
#                        selectInput("min_gtex_tpm", 
#                                    label = "Minimum TPM of non-placental tissues:", 
#                                    choices=list(`>0.1`=0.1,`>1`=1,`>5`=5, `>10`=10,`>100`=100),
#                                    selected=1),
#                        # drop down - min TPM of placenta 
#                        selectInput("min_gtex_fc", 
#                                    label = "Minimum fold change of a non-placental tissue compared with the placenta (i.e. TPM (non-placenta) / TPM (Placenta)):", 
#                                    choices=list(`>2x`=2,`>3x`=3,`>5x`=5,`>10x`=10,`>50x`=50,`>100x`=100),
#                                    selected=2),
#                        downloadButton("download_not_in_pt", "Download")
#                    ), # end of sidebarPanel
#                    # Show a plot of the generated distribution
#                    mainPanel(
#                        tabsetPanel(
#                            tabPanel("Summary",value="summary",
#                                        DT::dataTableOutput('not_in_placenta_summary')
#                            ),
#                            tabPanel("Rank",value="rank",
#                                     DT::dataTableOutput('not_in_placenta_rank')
#                            ),
#                            tabPanel("GO annotation",value="go.annotation",
#                                     DT::dataTableOutput('not_in_placenta_go')
#                            ),
#                            id="not.in.pt.tab"
#                        )
#                    ) # end of mainPanel
#                ) # end of sidebarLayout
#            ) # end of fluidPage
#        ) # end of tabPanel - not in placenta 

    ), # end of navbarMenu (Placenta vs GTEx)


    tabPanel(title="Genome Browser",
        fluidPage(
            #tag$head(tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js")) # does not work
            #withTags({
            #    head(
            #            script(
            #                src="http://www.biodalliance.org/release-0.13/dalliance-compiled.js"
            #            )
            #    )
            #}), # it works!
            #tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js"), # it works locally, but not shinyapps.io
            #tags$script(src = "dalliance-compiled.js"), # (www/dalliance-compiled.js) it works both locally and shinyapps.io
            tags$script(src = "dalliance-all.js"), # github version (0.13.7-dev)

            includeScript(path = "js/dalliance_ui.js"),

            tags$div(id="svgHolder","Dalliance goes here...")
            #tags$script(HTML("if (window.innerHeight < 400) alert('Screen too small');"))
        )
    ),

    #tabPanel(title="Download",
    #    "To Be Made..."
    #),

    tabPanel(title="TrackHub",
        includeMarkdown("trackhub.md")
    ),

    tabPanel("About",
        fluidPage(
            fluidRow(
                column(8,
                    includeMarkdown("about.md")
                ),
                column(4,
                    a("Tweets by ObsGynae", class="twitter-timeline", "data-theme"="light", "data-link-color"="#19CF86", href="https://twitter.com/ObsGynaeCam?ref_src=twsrc%5Etfw"),
                    tags$script(src = "https://platform.twitter.com/widgets.js", charset="utf-8")
                )
            ) # end of fluidRow
        )
    )
)
