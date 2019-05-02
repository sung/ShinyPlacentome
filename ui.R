# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019
# Last modified 1/May/2019

library(shiny)
#library(shinythemes)
# Define UI for application that draws a histogram
navbarPage(
    title="POPS Placenta Transcriptome",

    navbarMenu("DEG",
        tabPanel("By gene name",
            fluidPage(
                sidebarLayout(
                    # drop down 
                    sidebarPanel(
                        textInput("gene", "Gene Name (HGNC):", "FSTL3")
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
                    ) 
                )
            )
        ),
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
        )# end of tabPanel by FC
    ), # end of navbarMenu DEG
    
    navbarMenu("Search",
        tabPanel("ENSEMBL ID",
            fluidPage(
                sidebarLayout(
                    # drop down 
                    sidebarPanel(
                        textInput("ensg_id", "Text input:", "general")
                    ),
            
                    # Show a plot of the generated distribution
                    mainPanel("hello")
                )
            )
        )
    ),
    tabPanel("About",
            "This panel is intentionally left blank"
    )
)
