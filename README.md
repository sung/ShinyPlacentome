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
[1] nvimcom_0.9-83    data.table_1.12.2 colorout_1.2-2   

loaded via a namespace (and not attached):
[1] compiler_3.6.1 tools_3.6.1   
```

## How to start locally
```R
# run this command from your R session
shiny::runGitHub("ShinyPlacentome", "sung") 
```

----
Developed by Sung Gong <ssg29@cam.ac.uk>
