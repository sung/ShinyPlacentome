# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019

library(shiny)
library(DT)
library(data.table)

load("RData/DEG.RData")

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
    
    output$deg_pval_PE <- DT::renderDataTable(
        DT::datatable(dt.deseq[tr.type==input$transcript_pval & new.padj.x<input$pval,
                               .(Gene=hgnc_symbol,ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`log2FC(DESeq2)`=round(log2FoldChange.x,2),
                                 `P-val(DESeq2)`=round(pvalue.x,3),`adjusted pval`=round(new.padj.x,3))], 
                      rownames = FALSE, filter='top', options = list(pageLength = 15))
    )

    output$deg_pval_SGA <- DT::renderDataTable(
        DT::datatable(dt.deseq[tr.type==input$transcript_pval & new.padj.y<input$pval,
                               .(Gene=hgnc_symbol,ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`log2FC(DESeq2)`=round(log2FoldChange.y,2),
                                 `P-val(DESeq2)`=round(pvalue.y,3),`adjusted pval`=round(new.padj.y,3))], 
                      rownames = FALSE, filter='top', options = list(pageLength = 15))
    )
    
    output$deg_boot_oneThird<- DT::renderDataTable(
        DT::datatable(dt.boot[tr.type==input$transcript_boot & analysis.type=="oneThird",
                              .(Gene=hgnc_symbol,ID,"log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))],
                      rownames = FALSE, filter='top', options = list(pageLength = 15))
    )
    
    output$deg_boot_all<- DT::renderDataTable(
        DT::datatable(dt.boot[tr.type=="total-RNA:protein_coding" & analysis.type=="all",
                              .(Gene=hgnc_symbol,ID,"log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))],
                      rownames = FALSE, filter='top', options = list(pageLength = 15))
    )
    
    # By Gene Name(s) 
    gene.names<-dt.deseq[!is.na(hgnc_symbol),.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol
    updateSelectizeInput(session, 'genes', choices = gene.names, server = TRUE)

    output$deg_gene_pval <- DT::renderDataTable(
        DT::datatable(
            rbind(
                dt.deseq[hgnc_symbol %in% input$genes,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`adjusted pval`=round(new.padj.x,3))], 
                dt.deseq[hgnc_symbol %in% input$genes,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`adjusted pval`=round(new.padj.y,3))]
                )
        )
    )
    output$deg_gene_boot <- DT::renderDataTable(
        DT::datatable(
                dt.boot[hgnc_symbol %in% input$genes,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    )

    # By ENSG ID (s) 
    ensg.ids<-dt.deseq[grepl("^ENSG",ID),.N,ID][order(ID)]$ID#
    updateSelectizeInput(session, 'ensgs', choices = ensg.ids, server = TRUE)

    output$deg_ensg_pval <- DT::renderDataTable(
        DT::datatable(
            rbind(
                dt.deseq[ID %in% input$ensgs,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`adjusted pval`=round(new.padj.x,3))], 
                dt.deseq[ID %in% input$ensgs,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`adjusted pval`=round(new.padj.y,3))]
                )
        )
    )
    output$deg_ensg_boot <- DT::renderDataTable(
        DT::datatable(
                dt.boot[ID %in% input$ensgs,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    )
})
