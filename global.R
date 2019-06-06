gtex_tissues<-list(
`Adipose Tissue`="Adipose_Tissue",
`Adrenal Gland`="Adrenal_Gland",
`Blood`="Blood",
`Blood Vessel`="Blood_Vessel",
`Brain`="Brain",
`Breast`="Breast",
`Colon`="Colon",
`Esophagus`="Esophagus",
`Heart`="Heart",
`Liver`="Liver",
`Lung`="Lung",
`Muscle`="Muscle",
`Nerve`="Nerve",
`Pancreas`="Pancreas",
`Pituitary`="Pituitary",
`Small Intestine`="Small_Intestine",
`Skin`="Skin",
`Stomac`="Stomach",
`Spleen`="Spleen",
`Thyroid`="Thyroid")

# Suppl information from: https://academic.oup.com/bib/article/18/2/205/2562739
# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues
fTau <- function(x){
    if(all(!is.na(x))){
        if(min(x, na.rm=TRUE) >= 0){
            if(max(x)!=0){
                x <- (1-(x/max(x)))
                res <- sum(x, na.rm=TRUE)
                res <- res/(length(x)-1)
            }else{
                res <- 0
            }
 		}else{
            res <- NA
            #print("Expression values have to be positive!")
 		} 
 	}else{
        res <- NA
        #print("No data for this gene avalable.")
    } 
    return(res)
}

theme_Publication <- function(base_size=14, base_family="") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",size = rel(1.3), hjust = 0.5),
			text = element_text(),
			panel.background = element_rect(colour = NA),
			plot.background = element_rect(colour = NA),
			#panel.border = element_rect(colour = NA),
			axis.title = element_text(face = "bold",size = rel(1.2)),
			axis.title.y = element_text(angle=90,vjust =2),
			axis.title.x = element_text(vjust = -0.2),
			axis.text = element_text(size=rel(1.1)), 
			axis.line = element_line(colour="black"),
			axis.ticks = element_line(),
			panel.grid.major = element_line(colour="#f0f0f0"),
			panel.grid.minor = element_blank(),
			legend.key = element_rect(colour = NA),
			#legend.background = element_rect(colour = 'black', fill='grey',linetype='dashed'),
			legend.position = "right",
			#legend.direction = "horizontal",
			legend.key.size= unit(0.7, "cm"),
			legend.spacing= unit(0.1, "mm"),
			legend.title = element_text(face="bold.italic",size=rel(1)),
			legend.text = element_text(size = rel(0.9),family = "sans"),
			plot.margin=unit(c(10,5,5,5),"mm"),
			strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
			strip.text = element_text(face="bold",size=rel(1.1))
          ) + if(packageVersion("ggplot2")<=2.1){theme(legend.margin = unit(0.1, "mm"))}else{theme(legend.spacing = unit(0.1, "mm"))}
	   )
}

#enableBookmarking(store = "url")
enableBookmarking(store = "server")
