# <img width="410" height="100" src="https://user-images.githubusercontent.com/25957174/68218912-0faa6480-0020-11ea-9284-84f2bd761d4f.png"/>

# Geng_2019  
This repository contains the code for analysis of Hi-C data of the paper:
## "Three-dimensional Genome Structure Reveals Distinct Chromatin Signatures of Mouse Female Germline Stem Cells during Development"
Geng G. Tian, Xinyan Zhao, Wenhai Xie, Xiaoyong Li, Changliang Hou, Yinjuan Wang, Lijuan Wang, Xiaodong Zhao, Hua Li, Jing Li, and Ji Wu
### Installation
#### 1. Requirement  
* HiC-Pro:https://github.com/nservant/HiC-Pro  
* PSYCHIC:https://github.com/dhkron/PSYCHIC  
* R>3.4.0  
* perl  
* python  
* HiTC  
* Unix tools - cut, sed, pushd, popd (typically installed by default)  
* matrix2insulation.pl: https://github.com/dekkerlab/crane-nature-2015  
* samtools (>1.1）  
* bedtools

### Input matrix  
input_matrix should be in tab-delimitered file, which are came from the output of HiC-Pro. The contact maps are then available in the hic_results/matrix folder. Raw contact maps are in the raw folder and normalized contact maps in the iced folder. A contact map is defined by : A list of genomic intervals related to the specified resolution (BED format). A matrix, stored as standard triplet sparse format (i.e. list format).

### Session Info
```R
> sessionInfo()
R version 3.5.2 (2018-12-20)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.2 LTS

Matrix products: default
BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
[1] ggplot2_3.2.1        reshape2_1.4.3       limma_3.38.3
[4] HiTC_1.26.0          GenomicRanges_1.34.0 GenomeInfoDb_1.18.2
[7] IRanges_2.16.0       S4Vectors_0.20.1     BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0                  compiler_3.5.2
 [3] pillar_1.3.1                RColorBrewer_1.1-2
 [5] plyr_1.8.4                  XVector_0.22.0
 [7] bitops_1.0-6                tools_3.5.2
 [9] zlibbioc_1.28.0             tibble_2.0.1
[11] gtable_0.2.0                lattice_0.20-38
[13] pkgconfig_2.0.2             rlang_0.4.0
[15] Matrix_1.2-17               DelayedArray_0.8.0
[17] GenomeInfoDbData_1.2.0      withr_2.1.2
[19] dplyr_0.8.0.1               rtracklayer_1.42.1
[21] stringr_1.3.1               Biostrings_2.50.2
[23] tidyselect_0.2.5            grid_3.5.2
[25] glue_1.3.0                  Biobase_2.42.0
[27] R6_2.3.0                    XML_3.98-1.17
[29] BiocParallel_1.16.6         purrr_0.3.0
[31] magrittr_1.5                Rsamtools_1.34.1
[33] scales_1.0.0                matrixStats_0.54.0
[35] GenomicAlignments_1.18.1    assertthat_0.2.0
[37] SummarizedExperiment_1.12.0 colorspace_1.4-0
[39] stringi_1.2.4               RCurl_1.95-4.11
[41] lazyeval_0.2.1              munsell_0.5.0
[43] crayon_1.3.4
 ```
### Bugs and Feedback
#### For bugs, questions and discussions please use the Github Issues or send me the E-mail: gengtian@sjtu.edu.cn.
