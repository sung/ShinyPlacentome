# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019
# Last modified 1/May/2019

library(shiny)
#library(shinythemes)
# Define UI for application that draws a histogram

#load("RData/DEG.RData") # client-side gene selection with selectize
#gene.names<-dt.deseq[!is.na(hgnc_symbol),.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol # client side list of genes

navbarPage(
#    withTags({
#        head(
#                script(
#                    src="http://www.biodalliance.org/release-0.13/dalliance-compiled.js"
#                )
#        )
#    }),

    title="POPS Placenta Transcriptome",

    #tag$head(tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js")), # does not work

    navbarMenu("Browse DEG",
        tabPanel("By P-value",
            fluidPage(
                sidebarLayout(
                    sidebarPanel(
                        # drop down 
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
                            tabPanel("One third of qualified genes",
                                     DT::dataTableOutput('deg_boot_oneThird')
                            ),
                            tabPanel("All qualified genes",
                                     DT::dataTableOutput('deg_boot_all')
                            )
                        )
                    ) # end of mainPanel
                ) # end of sidebarLayout
            ) # end of fluidPage
        ) # end of tabPanel by FC
    ), # end of navbarMenu DEG
    
    navbarMenu("Search",
        tabPanel("By gene name",
            fluidPage(
                sidebarLayout(
                    # drop down 
                    sidebarPanel(
                    #    textInput("gene", "Gene Name (HGNC):", "FSTL3") # old-school
                    #    selectizeInput("genes","Gene Name(s):", choices=gene.names, selected = "FSTL3", multiple=TRUE) # client side
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
        ),
        tabPanel("By ENSEMBL ID",
            fluidPage(
                sidebarLayout(
                    # drop down 
                    sidebarPanel(
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
        )
    ),

    tabPanel(title="Genome Browser",
        #tag$head(tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js")) # does not work
        tags$script(src = "http://www.biodalliance.org/release-0.13/dalliance-compiled.js") # it works
        #withTags({
        #    head(
        #            script(
        #                src="http://www.biodalliance.org/release-0.13/dalliance-compiled.js"
        #            )
        #    )
        #}) # this also works!
    ),

    tabPanel("About",
            "This panel is intentionally left blank"
    )
)
