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
    # below works
#    withTags({
#        head(
#                script(
#                    src="http://www.biodalliance.org/release-0.13/dalliance-compiled.js"
#                )
#        )
#    }),

    #tag$head(tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js")), # does not work

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
                                     DT::dataTableOutput('deg_boot_pe_all'),
                                     hr(),
                                     DT::dataTableOutput('deg_boot_pe_oneThird')
                            ),
                            tabPanel("SGA",
                                     DT::dataTableOutput('deg_boot_sga_all'),
                                     hr(),
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
                    #    selectizeInput("genes","Gene Name(s):", choices=gene.names, selected = "FSTL4", multiple=TRUE) # client side
                        selectizeInput("genes","Gene Name(s):", choices=NULL, selected="FSTL3", multiple=TRUE) # server side
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
                        selectizeInput("ensgs","ENSG ID(s):", choices=NULL, multiple=TRUE) # server side
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
    
    tabPanel(title="Abundance Level",
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
                                        "lincRNA"="lincRNA",                # total-RNA
                                        "pseudogene"="pseudogene",          # total-RNA
                                        "circRNA"="circRNA",
                                        "miRNA"="miRNA",
                                        "piRNA"="piRNA",
                                        "novel transcript isoforms"="novel_isoform",
                                        "novel miRNA"="novel_miRNA",
                                        "novel small-RNA"="novel_smallRNA"),
                            selected="protein_coding"),
                    # conditional drop down - min FPKM of placenta (except circRNA)
                    conditionalPanel(
                        condition="input.ab_transcript != 'circRNA'",
                        selectInput("fpkm", 
                                    label = "Minimum FPKM:", 
                                    choices=list(`>0.01`=0.01,`>0.1`=0.1,`>1`=1,`>5`=5,`>10`=10),
                                    selected=1)),
                    # conditional radio button - polyA+ or not (circRNA only)
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
                    # a set of radio buttons - transcript type
                ),
                # Show a plot of the generated distribution
                mainPanel(
                    DT::dataTableOutput('pops_tr')
                ) # end of mainPanel
            ) # end of sidebarLayout
        ) # end of fludPage
    ), # end of tabPanel - reconstructed transcriptome

    tabPanel(title="Placenta specific",
        fluidPage(
            sidebarLayout(
                sidebarPanel(
                    helpText("Browse genes expressed specifically in the placenta"),
                    # drop down 
                    selectInput("transcript_tau", 
                            label="Choose a type of transcript:",
                            choices = list(
                                        "protein coding"="protein_coding", 
                                        "lincRNA"="lincRNA",
                                        "processed pseudogene"="processed_pseudogene")),
                    # drop down - min FPKM of placenta 
                    selectInput("pt_fpkm", 
                                label = "Minimum FPKM of Placenta:", 
                                choices=list(`>0.1`=0.1,`>1`=1,`>5`=5, `>10`=10),
                                selected=1),
                    # a set of radio buttons - transcript type
                    # tau score range 
                    sliderInput("tau", label = "Tau score", min = 0.9, max = 1, value = c(0.99,1)),
                    # drop down - min FPKM of placenta 
                    selectInput("pt_gtex_fc", 
                                label = "Fold change of placenta compared with the average of 20 GTEx tissues:", 
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
                    d3heatmapOutput("heatmap", width="90%", height="1200px"),
                    hr(),
                    DT::dataTableOutput('tau')
                ) # end of mainPanel
            ) # end of sidebarLayout
        ) # end of fludPage
    ), # end of tabPanel

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
            tags$script(src = "dalliance-compiled.js"), # (www/dalliance-compiled.js) it works both locally and shinyapps.io

            includeScript(path = "js/dalliance_ui.js"),

            tags$div(id="svgHolder","Dalliance goes here...")
            #tags$script(HTML("if (window.innerHeight < 400) alert('Screen too small');"))
        )
    ),

    #tabPanel(title="Download",
    #    "To Be Made..."
    #),

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
