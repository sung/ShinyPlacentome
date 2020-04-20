# The landscape of placenta transcriptome in health and disease
# Sung Gong <ssg29@cam.ac.uk>
# https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/
# https://github.com/sung
# First created 1/May/2019

library(data.table)
library(DT)
library(d3heatmap)
library(ggplot2)
library(ggsci)
library(shiny)

load("RData/DEG.RData") # dt.deseq (isa data.table) # dt.boot (isa data.table) 3.1M
load("RData/dt.gtex.pt.tpm.tau.RData") # dt.gtex.pt.tpm.tau (isa data.table) 3.2M # GRCh38.90
#load("RData/dt.pops.tr.RData") # dt.pops.tr # it takes long
#library(feather) # devtools::install_github("wesm/feather/R") 
                 # https://blog.rstudio.com/2016/03/29/feather/
#dt.pops.tr = data.table(feather::read_feather("RData/dt.pops.tr.feather")) # file too big (473MB)
load("RData/dl.abundance.RData") # bin/R/Placentome/dl.pops.tr.abundance.R # 5.2M
load("RData/dt.gtex.tpm.RData") # 6.2M (ensembl_gene_id,hgnc_symbol,gene_biotype,baseMean,TPM,Tissue) # GRCh38.90
load("RData/dt.ensg.desc.2019-05-20.RData") # 1.3M
load("RData/dt.ensg.go.2019-05-22.RData") # 5.1M

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
    ##############
    # Browse DEG #
    ##############
    # 1. by p-value
    dt.deseq.pe<-reactive({
        dt.deseq[tr.type==input$transcript_pval & new.padj.x<as.numeric(input$pval),.(Gene=hgnc_symbol,ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))]

    })

    dt.deseq.sga<-reactive({
        dt.deseq[tr.type==input$transcript_pval & new.padj.y<as.numeric(input$pval),.(Gene=hgnc_symbol,ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
    })

    output$deg_pval_PE <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.deseq.pe(), rownames = FALSE, filter='top', options = list(pageLength = 15))
    })
    output$deg_pval_SGA <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
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
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.boot.pe()[analysis.type=="all"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from all the qualified genes", rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$deg_boot_pe_oneThird <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.boot.pe()[analysis.type=="oneThird"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from one third of highly abudant genes",rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$deg_boot_sga_all <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.boot.sga()[analysis.type=="all"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from all the qualified genes", rownames = FALSE, filter='top', options = list(pageLength = 15))
    })
    output$deg_boot_sga_oneThird <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.boot.sga()[analysis.type=="oneThird"][,-"analysis.type"], caption="Table 1. Top 5% (fold-chage) selected from one third of highly abudant genes",rownames = FALSE, filter='top', options = list(pageLength = 15))
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

    
    # 3. DEG - By Gene Name(s) 
    deg.gene.names<-dt.deseq[!is.na(hgnc_symbol),.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol
    updateSelectizeInput(session, 'deg_genes', choices = deg.gene.names, server = TRUE)

    output$deg_gene_pval <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        DT::datatable(
            rbind(
                dt.deseq[hgnc_symbol %in% input$deg_genes,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))], 
                dt.deseq[hgnc_symbol %in% input$deg_genes,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
                )
        )
    })
    output$deg_gene_boot <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(
                dt.boot[hgnc_symbol %in% input$deg_genes,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"Log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"Log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    })

    # 4. DEG - By ENSG ID (s) 
    ensg.ids<-dt.deseq[grepl("^ENSG",ID),.N,ID][order(ID)]$ID#
    updateSelectizeInput(session, 'deg_ensgs', choices = ensg.ids, server = TRUE)

    output$deg_ensg_pval <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(
            rbind(
                dt.deseq[ID %in% input$deg_ensgs,
                         .(`Subject`="PE",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.x,1),FPKM=round(meanFpkm.x,2),`Log2FC(DESeq2)`=round(log2FoldChange.x,2),`P-val(DESeq2)`=round(pvalue.x,3),`Adjusted P-val`=round(new.padj.x,3))], 
                dt.deseq[ID %in% input$deg_ensgs,
                         .(`Subject`="SGA",`Gene`=hgnc_symbol, `ENSG`=ID,Count=round(baseMean.y,1),FPKM=round(meanFpkm.y,2),`Log2FC(DESeq2)`=round(log2FoldChange.y,2),`P-val(DESeq2)`=round(pvalue.y,3),`Adjusted P-val`=round(new.padj.y,3))]
                )
        )
    })
    output$deg_ensg_boot <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(
                dt.boot[ID %in% input$deg_ensgs,
                        .(`Analysis type`=analysis.type,`Gene`=hgnc_symbol,`ENSG`=ID,"Log2FC(PE vs Control)"=round(log2FC.Boot.x,2),"Log2FC(SGA vs Control)"=round(log2FC.Boot.y,2))]
        )
    })


    ###############
    ## Abundance ##
    ###############
    # 1. prepare data.table 
    dt.abundance<-reactive({
        if(input$ab_transcript=="circRNA"){
            if(input$in_polyA){
                dl.abundance[[input$ab_transcript]][!(`In PolyA+?`) & freq>input$evi.ratio[1] & freq<=input$evi.ratio[2]]
            }else{
                dl.abundance[[input$ab_transcript]][freq>input$evi.ratio[1] & freq<=input$evi.ratio[2]]
            }
        }else if(input$ab_transcript=="novel_isoform"){
            num_exon<-ifelse(input$no_single_exon,1,0)
            if(as.numeric(input$fpkm)==0){
                dl.abundance[[input$ab_transcript]][exon.cnt>num_exon & `sample freq`>=input$evi.ratio[1] & `sample freq`<=input$evi.ratio[2],-"class_code"]
            }else{
                dl.abundance[[input$ab_transcript]][exon.cnt>num_exon & FPKM>as.numeric(input$fpkm) & `sample freq`>=input$evi.ratio[1] & `sample freq`<=input$evi.ratio[2],-"class_code"]
            }
        }else{
            if(as.numeric(input$fpkm)==0){
                dl.abundance[[input$ab_transcript]]
            }else{
                dl.abundance[[input$ab_transcript]][FPKM>as.numeric(input$fpkm)]
            }
        }
    })

    # 2. render data.table
    output$pops_tr <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.abundance(),rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    # 3. download bits
    output$download_pops_tr<- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste("POPS.transcriptome",input$ab_transcript,"csv.gz", sep = ".")
        },

        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            write.csv(dt.abundance(), gzfile(file), row.names = FALSE, quote=F)
        }
    ) # end of downloadData

    # 4. ggbio for the Manhattan plot
    gr.abundance<-reactive({
        dt.manhattan<-dt.abundance()
        dt.manhattan[,`score`:=ifelse(Type=="circRNA",meanCount/10^3,meanFpkm/10^6)]
        GenomicRanges::makeGRangesFromDataFrame(dt.manhattan,keep.extra.columns=T)
    })

    #######################
    ## Placenta-specific ##
    #######################
    #output$options<- renderPrint({ paste(input$transcript_tau, input$pt_tpm1, input$tau[1], input$tau[2], input$pt_gtex_fc)})
    #output$test4<- renderPrint({ lapply(input, class)})
    # 1. tau-based
    dt.tau<-reactive({
        dt.gtex.pt.tpm.tau[!grepl("^HIST",hgnc_symbol) & !hgnc_symbol%in%c("RMRP","AL356488.2") &
                            gene_biotype==input$transcript_tau &  Placenta>as.numeric(input$pt_tpm1) & Tau>input$tau[1] & Tau<=input$tau[2] & Placenta/meanTPMGTEx > as.numeric(input$pt_gtex_fc),-"meanTPMGTEx"]
    })

    row.num<-reactive({
        ifelse(nrow(dt.tau())>50,50,nrow(dt.tau()))
    })

    mat.tau<-reactive({
        dt.tau<-dt.tau()
        mat.tau<-as.matrix(dt.tau[,-c("ensembl_gene_id","Tau","chromosome_name","hgnc_symbol","gene_biotype")])[1:row.num(),]
        mat.tau<-log10(mat.tau+0.001)
        colnames(mat.tau)<-gsub("_"," ",colnames(mat.tau))
        rownames(mat.tau)<-dt.tau[1:row.num()]$hgnc_symbol
        return(mat.tau)
    })

    output$tau <- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.tau()[,-"gene_biotype"], caption="Table 1. Expression level (TPM) of 20 somatic tissues (from GTEx) and the placenta (this study)", rownames = FALSE, filter='top', options = list(pageLength = 15))
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
        paste("Clustering of top ", row.num(), input$transcript_tau, "genes specifically expressed in the placneta (color-scale:log10(TPM))")
    })

    output$heatmap <- renderD3heatmap({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Rendering heatmap", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        d3heatmap(
            mat.tau(),
            cellnote=dt.tau()[,-c("ensembl_gene_id","Tau","chromosome_name","hgnc_symbol","gene_biotype")][1:row.num()],
            dendrogram = "column",
            xaxis_font_size = "10pt",
            colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name="RdYlBu")))(100) 
        )
    })

    # 2. ALL or by Gene Name(s) 
    gtex.gene.names<-dt.gtex.pt.tpm.tau[,.N,hgnc_symbol][order(hgnc_symbol)]$hgnc_symbol
    updateSelectizeInput(session, 'gtex_genes', choices = gtex.gene.names, server = TRUE)

    dt.gtex<-reactive({
        if(input$radio_gtex==1){ # all the genes
            dt.gtex.tau<-cbind(dt.gtex.pt.tpm.tau[,1:4], dt.gtex.pt.tpm.tau[,.(Tau)],dt.gtex.pt.tpm.tau[,5:25])
            dt.gtex.tau[Placenta >= as.numeric(input$pt_tpm2)][order(-Placenta)]
        }else{ # user-provided genes
            dt.foo<-melt.data.table(dt.gtex.pt.tpm.tau[hgnc_symbol %in% input$gtex_genes, -c("meanTPMGTEx","Tau")],
                            id.vars=c("chromosome_name","ensembl_gene_id","hgnc_symbol","gene_biotype"), variable.name="Tissue", value.name="TPM")
            dt.foo[,`Source`:=ifelse(Tissue=="Placenta","POPS","GTEx")][order(hgnc_symbol,-Source,-TPM)]
        }
    })

    output$gtex_tpm<- DT::renderDataTable({
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        DT::datatable(dt.gtex(), caption="Table 1. Expression level (TPM) of 20 somatic tissues (from GTEx) and the placenta (this study)", rownames = FALSE, filter='top', options = list(pageLength = 21))
    })

    output$gtex_tpm_barchart<-renderPlot({
        if(length(input$gtex_genes>0)){
            ggplot(dt.gtex(), aes(Tissue, log10(TPM+1), fill=hgnc_symbol)) + 
                geom_bar(col='gray10',stat="identity",position="dodge") + 
                scale_x_discrete(limits=dt.gtex()[,.(`meanTPM`=mean(TPM)),Tissue][order(-`meanTPM`)]$Tissue) +
                ggsci::scale_fill_jco(name="Gene Name(s)",alpha=.75) +
                ylab("log10(TPM+1)") +
                theme_Publication() + 
                theme(axis.text.x=element_text(angle=45, hjust=1))
        }
    })

    output$download_gtex<- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            if(input$radio_gtex==1){
                paste("TPM.GTEx.vs.placenta.more.than",input$pt_tpm2,"TPM.csv.gz", sep = ".")
            }else{
                "TPM.GTEx.vs.placenta.csv"
            }
        },
        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            if(input$radio_gtex==1){
                write.csv(dt.gtex(), gzfile(file), row.names = FALSE, quote=F)
            }else{
                write.csv(dt.gtex(), file, row.names = FALSE, quote=F)
            }
        }
    ) # end of downloadData

    #####################
    ## Not in Placenta ##
    #####################
    # 20 GTEx tissue + 'Placenta'
    # By deault, n=2('Breast' AND 'Blood') omiited from 20 GTEx
    dt.pt.bottom<-reactive({
        dt.gtex.desc<-merge(dt.gtex.tpm,dt.ensg.desc[,.(ensembl_gene_id,chromosome_name,description)],all.x=T)
        my.ensg<-dt.gtex.desc[!Tissue %in% c("Placenta",input$no_gtex_tissue)
                              & baseMean > as.numeric(input$min_gtex_count)
                              & !chromosome_name %in% c("Y")
                              & gene_biotype==input$transcript_not_in_pt
                              ,.N,ensembl_gene_id][N==length(gtex_tissues)-length(input$no_gtex_tissue),.N,ensembl_gene_id]$ensembl_gene_id
        dt.gtex.rank<-dt.gtex.desc[!Tissue %in% input$no_gtex_tissue 
                                   & ensembl_gene_id %in% my.ensg
                                   ,.(Tissue,TPM,rank=length(gtex_tissues)+1-length(input$no_gtex_tissue)+1 -rank(TPM))
                                   ,.(chromosome_name,ensembl_gene_id,hgnc_symbol,description)][order(ensembl_gene_id,rank)] 
        if(input$no_ribosomal){
            dt.bottom<-dt.gtex.rank[!grepl("ribosom",description) & rank==length(gtex_tissues)+1-length(input$no_gtex_tissue) & Tissue=="Placenta"][order(-TPM)] #
        }else{
            dt.bottom<-dt.gtex.rank[rank==length(gtex_tissues)+1-length(input$no_gtex_tissue) & Tissue=="Placenta"][order(-TPM)] #
        }
        # genes where the placenta ranked at the bottom
        dt.gtex.rank[ensembl_gene_id %in% dt.bottom$ensembl_gene_id][order(ensembl_gene_id,rank)]
    })

    # apply min TPM of GTEx & min FC
    dt.pt.bottom.summary<-reactive({
        dt.pt.bottom<-dt.pt.bottom()
        # get min,max,mean,median of the genes above
        dt.pt.bottom.summary<-merge(
                                    merge(
                                        dt.pt.bottom[,.(Tau=round(sapply(.SD,fTau),3)),.(chromosome_name,ensembl_gene_id,hgnc_symbol,description),.SDcol="TPM"],
                                        dt.pt.bottom[Tissue!="Placenta",
                                                    .(`GTEx_minTPM`=min(TPM),`GTEx_maxTPM`=max(TPM),`GTEx_medianTPM`=median(TPM),`GTEx_meanTPM`=round(mean(TPM),3)),
                                                    .(ensembl_gene_id)]
                                    ),
                                    dt.pt.bottom[Tissue=="Placenta",.(ensembl_gene_id,`Placenta_meanTPM`=TPM)]
                                          )[order(Placenta_meanTPM)]
        dt.pt.bottom.summary[GTEx_minTPM>as.numeric(input$min_gtex_tpm) & GTEx_minTPM/Placenta_meanTPM>as.numeric(input$min_gtex_fc)]
    })

    dt.pt.bottom.rank<-reactive({
        dt.pt.bottom()[ensembl_gene_id %in% dt.pt.bottom.summary()$ensembl_gene_id][order(ensembl_gene_id,rank)]
    })

    dt.pt.bottom.go<-reactive({
        if(TRUE){
            dt.ensg.go[ensembl_gene_id %in% dt.pt.bottom.summary()$ensembl_gene_id]
        }else{
            data.table(biomaRt::getBM(attributes =c("ensembl_gene_id","hgnc_symbol","description","go_id","name_1006"), 
                                                filters = "ensembl_gene_id", 
                                                values = dt.pt.bottom.summary()$ensembl_gene_id,
                                                mart=biomaRt::useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
                                                ))
        }
    })

    output$not_in_placenta_summary<- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        DT::datatable(dt.pt.bottom.summary(), rownames = FALSE, filter='top', options = list(pageLength = 15))
    })

    output$not_in_placenta_rank<- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        DT::datatable(dt.pt.bottom.rank()[,-"description"], rownames = FALSE, filter='top', options = list(pageLength = length(gtex_tissues)+1-length(input$no_gtex_tissue)))
    })

    output$not_in_placenta_go<- DT::renderDataTable({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Loading table", value = 0)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        DT::datatable(dt.pt.bottom.go(), rownames = FALSE, filter='top', 
                      options = list(searchHighlight = TRUE, search = list(regex = FALSE, caseInsensitive = TRUE, search = 'DNA repair'), pageLength = 15))
    })

    output$download_not_in_pt<- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
                paste(input$not.in.pt.tab,"of",input$transcript_not_in_pt,"genes.not.in.placneta.but.in.gtex.at.least",
                      input$min_gtex_count,"count",input$min_gtex_tpm,"tpm",input$min_gtex_fc,"FC.tsv.gz",sep=".")
        },
        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            if(input$not.in.pt.tab=="summary"){
                write.table(dt.pt.bottom.summary(),gzfile(file), sep="\t",row.names = FALSE, quote=F)
            }else if(input$not.in.pt.tab=="rank"){
                write.table(dt.pt.bottom.rank(),gzfile(file), sep="\t",row.names = FALSE, quote=F)
            }else{
                write.table(dt.pt.bottom.go(),gzfile(file), sep="\t",row.names = FALSE, quote=F)
            }
        }
    ) # end of downloadData

})
