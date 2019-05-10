# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019

library(data.table)
library(DT)
library(d3heatmap)
library(shiny)

load("RData/DEG.RData") # dt.deseq (isa data.table)
                        # dt.boot (isa data.table)
load("RData/dt.gtex.pt.fpkm.tau.RData") # dt.gtex.pt.fpkm.tau (isa data.table)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
    
    ##########
    # Browse #
    ##########
    # 1. by p-value
    dt.deseq.pe<-reactive({
        dt.deseq[tr.type==input$transcript_pval & new.padj.x<as.numeric(input$pval),.(Gene=hgnc_symbol,ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))]

    })

    dt.deseq.sga<-reactive({
        dt.deseq[tr.type==input$transcript_pval & new.padj.y<as.numeric(input$pval),.(Gene=hgnc_symbol,ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
    })

    output$deg_pval_PE <- DT::renderDataTable({
        DT::datatable(dt.deseq.pe(), rownames = FALSE, filter='top', options = list(pageLength = 15))
    })
    output$deg_pval_SGA <- DT::renderDataTable({
        DT::datatable(dt.deseq.sga(), rownames = FALSE, filter='top', options = list(pageLength = 15))
    })
    # downloadHandler() takes two arguments, both functions.
    # The content function is passed a filename as an argument, and
    # it should write out data to that filename.
    output$download_pval <- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste(input$transcript_pval, "p-value_less_than", input$pval,"csv", sep = ".")
        },

        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            write.csv(rbind(dt.deseq.pe()[,`Subject`:="PE"], dt.deseq.sga()[,`Subject`:="SGA"]), file, row.names = FALSE, quote=F)
        }
    ) # end of downloadData

    # 2. by FC (bootstrap) 
    dt.boot.pe<-reactive({
        dt.boot[tr.type==input$transcript_boot & !is.na(log2FC.Boot.x), .(analysis.type, Gene=hgnc_symbol,ID,"Log2FC(Boot)"=round(log2FC.Boot.x,3))]
    })

    dt.boot.sga<-reactive({
        dt.boot[tr.type==input$transcript_boot & !is.na(log2FC.Boot.y), .(analysis.type, Gene=hgnc_symbol,ID,"Log2FC(Boot)"=round(log2FC.Boot.y,3))]
    })

    output$deg_boot_pe_all <- DT::renderDataTable({
        DT::datatable(dt.boot.pe()[analysis.type=="all"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from all the qualified genes", rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$deg_boot_pe_oneThird <- DT::renderDataTable({
        DT::datatable(dt.boot.pe()[analysis.type=="oneThird"][,-"analysis.type"], caption="Table 2. Top 5% (fold-chage) selected from one third of highly abudant genes",rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$deg_boot_sga_all <- DT::renderDataTable({
        DT::datatable(dt.boot.sga()[analysis.type=="all"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from all the qualified genes", rownames = FALSE, filter='top', options = list(pageLength = 15))
    })
    output$deg_boot_sga_oneThird <- DT::renderDataTable({
        DT::datatable(dt.boot.sga()[analysis.type=="oneThird"][,-"analysis.type"], caption="Table 2. Top 5% (fold-chage) selected from one third of highly abudant genes",rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$download_boot<- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste(input$transcript_pval,"bootstrapping_analysis.csv", sep = ".")
        },

        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            write.csv(rbind(dt.boot.pe()[,`Subject`:="PE"], dt.boot.sga()[,`Subject`:="SGA"]), file, row.names = FALSE, quote=F)
        }
    ) # end of downloadData

    
    ##########
    # Search #
    ##########
    # By Gene Name(s) 
    gene.names<-dt.deseq[!is.na(hgnc_symbol),.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol
    updateSelectizeInput(session, 'genes', choices = gene.names, server = TRUE)

    output$deg_gene_pval <- DT::renderDataTable({
        DT::datatable(
            rbind(
                dt.deseq[hgnc_symbol %in% input$genes,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))], 
                dt.deseq[hgnc_symbol %in% input$genes,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
                )
        )
    })
    output$deg_gene_boot <- DT::renderDataTable({
        DT::datatable(
                dt.boot[hgnc_symbol %in% input$genes,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"Log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"Log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    })
    # By ENSG ID (s) 
    ensg.ids<-dt.deseq[grepl("^ENSG",ID),.N,ID][order(ID)]$ID#
    updateSelectizeInput(session, 'ensgs', choices = ensg.ids, server = TRUE)

    output$deg_ensg_pval <- DT::renderDataTable({
        DT::datatable(
            rbind(
                dt.deseq[ID %in% input$ensgs,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))], 
                dt.deseq[ID %in% input$ensgs,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
                )
        )
    })
    output$deg_ensg_boot <- DT::renderDataTable({
        DT::datatable(
                dt.boot[ID %in% input$ensgs,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"Log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"Log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    })

    #######################
    ## Placenta-specific ##
    #######################
    #output$options<- renderPrint({ paste(input$transcript_tau, input$pt_fpkm, input$tau[1], input$tau[2], input$pt_gtex_fc)})
    #output$test4<- renderPrint({ lapply(input, class)})

    dt.tau<-reactive({
            dt.gtex.pt.fpkm.tau[!grepl("^HIST",hgnc_symbol) & hgnc_symbol!="RMRP" &
                                gene_biotype==input$transcript_tau &  Placenta>as.numeric(input$pt_fpkm) & Tau>input$tau[1] & Tau<=input$tau[2] & Placenta/meanFpkmGTEx > as.numeric(input$pt_gtex_fc)]
    })

    row.num<-reactive({
        ifelse(nrow(dt.tau())>50,50,nrow(dt.tau()))
    })

    mat.tau<-reactive({
        dt.tau<-dt.tau()
        mat.tau<-as.matrix(dt.tau[,-c("ensembl_gene_id","meanFpkmGTEx","Tau","chromosome_name","hgnc_symbol","gene_biotype")])[1:row.num(),]
        mat.tau<-log10(mat.tau+0.001)
        colnames(mat.tau)<-gsub("_"," ",colnames(mat.tau))
        rownames(mat.tau)<-dt.tau[1:row.num()]$hgnc_symbol
        return(mat.tau)
    })

    output$tau <- DT::renderDataTable({
        DT::datatable(dt.tau()[,-"gene_biotype"], caption="Table 1. Expression level (FPKM) of 20 somatic tissues (from GTEX) and the placenta (this study)", rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$download_tau<- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste("placenta_specific",input$transcript_tau,"csv", sep = ".")
        },

        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            write.csv(dt.tau(), file, row.names = FALSE, quote=F)
        }
    ) # end of downloadData

    output$heatmap_title <- renderText({
        paste("Clustering of top ", row.num(), input$transcript_tau, "genes specifically expressed in the placneta (color-scale:log10(FPKM))")
    })

    output$heatmap <- renderD3heatmap({
        d3heatmap(
            mat.tau(),
            cellnote=dt.tau()[,-c("ensembl_gene_id","meanFpkmGTEx","Tau","chromosome_name","hgnc_symbol","gene_biotype")][1:row.num()],
            dendrogram = "column",
            xaxis_font_size = "10pt",
            colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name="RdYlBu")))(100) 
        )
    })


})
