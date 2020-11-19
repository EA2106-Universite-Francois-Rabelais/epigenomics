# Epigenomics

Scripts and data to reproduce results presented in this [article](https://www.mdpi.com/1422-0067/21/17/6028).
Differential methylation and differential expression are respectively detected with the DSS and edgeR R package. Data are processed with the data.table package for fast and efficient computations.

[](data/methylation_profiles.png)

## Required R packages

```
library(data.table)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(FactoMineR)
library(edgeR)
library(dynamicTreeCut)
library(pheatmap)
```

## Load data and R objects

## Data exploration

## Session info

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dynamicTreeCut_1.63-1 edgeR_3.26.8          limma_3.40.0          FactoMineR_2.3        ggdendro_0.1.22      
[6] reshape2_1.4.4        ggplot2_3.3.2         data.table_1.13.2    

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5           pillar_1.4.6         compiler_3.6.3       plyr_1.8.6           tools_3.6.3         
 [6] digest_0.6.27        lifecycle_0.2.0      tibble_3.0.4         gtable_0.3.0         lattice_0.20-41     
[11] pkgconfig_2.0.3      rlang_0.4.8          rstudioapi_0.11      ggrepel_0.8.2        yaml_2.2.1          
[16] withr_2.3.0          dplyr_1.0.2          stringr_1.4.0        cluster_2.1.0        generics_0.1.0      
[21] vctrs_0.3.4          locfit_1.5-9.4       flashClust_1.01-2    grid_3.6.3           tidyselect_1.1.0    
[26] scatterplot3d_0.3-41 glue_1.4.2           R6_2.5.0             farver_2.0.3         purrr_0.3.4         
[31] magrittr_1.5         scales_1.1.1         ellipsis_0.3.1       MASS_7.3-53          leaps_3.1           
[36] colorspace_1.4-1     labeling_0.4.2       stringi_1.5.3        munsell_0.5.0        crayon_1.3.4  
```
