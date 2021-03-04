[![DOI](https://zenodo.org/badge/184614898.svg)](https://zenodo.org/badge/latestdoi/184614898)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
## Introduction
A ShinyApp to search differentially expressed genes in the human placenta.
This is a source code repository for https://www.obgyn.cam.ac.uk/placentome/.

## Dependencies
1. DT
2. data.table
3. markdown
4. d3heatmap
5. RColorBrewer
6. grDevices

## `SessionInfo()` shown below (or [here](sessionInfo.txt))
```r
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] markdown_1.1      shinythemes_1.1.2 shiny_1.3.2       ggsci_2.9         ggplot2_3.1.1     d3heatmap_0.6.1.2
 [7] DT_0.6            nvimcom_0.9-83    data.table_1.12.2 colorout_1.2-2   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       pillar_1.4.0     compiler_3.6.1   later_1.0.0      plyr_1.8.4       base64enc_0.1-3 
 [7] tools_3.6.1      digest_0.6.22    tibble_2.1.1     gtable_0.3.0     pkgconfig_2.0.2  png_0.1-7       
[13] rlang_0.4.1      withr_2.1.2      dplyr_0.8.1      htmlwidgets_1.3  grid_3.6.1       tidyselect_0.2.5
[19] glue_1.3.1       R6_2.4.0         purrr_0.3.3      magrittr_1.5     scales_1.0.0     promises_1.1.0  
[25] htmltools_0.4.0  assertthat_0.2.1 mime_0.7         colorspace_1.4-1 xtable_1.8-4     httpuv_1.5.2    
[31] lazyeval_0.2.2   munsell_0.5.0    crayon_1.3.4    
```

## How to start locally
```R
# run this command from your R session
shiny::runGitHub("ShinyPlacentome", "sung") 
```

----
Developed by Sung Gong <ssg29@cam.ac.uk>
