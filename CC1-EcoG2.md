R Notebook
================

#Methods

##Amplicon bioinformatics: from raw reads to tables

``` r
.cran_packages <- c("ggplot2", "gridExtra", "devtools")
install.packages(.cran_packages) 
```

    ## Installing packages into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
BiocManager::install(.bioc_packages)
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## '?repositories' for details
    ## 
    ## replacement repositories:
    ##     CRAN: https://packagemanager.rstudio.com/all/__linux__/focal/latest

    ## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'dada2' 'phyloseq' 'DECIPHER' 'phangorn'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     Matrix

``` r
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ## Loading required package: ggplot2

    ## Loading required package: gridExtra

    ## Loading required package: devtools

    ## Loading required package: usethis

    ## Loading required package: dada2

    ## Loading required package: Rcpp

    ## Loading required package: phyloseq

    ## Loading required package: DECIPHER

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

    ## Loading required package: phangorn

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

    ##   ggplot2 gridExtra  devtools     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
set.seed(100)
miseq_path <- "/home/rstudio/MiSeq_SOP"
list.files(miseq_path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

##Filter and Trim

``` r
fnFs <- sort(list.files(miseq_path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern = "_R2_001.fastq"))
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#file.path permet d'adapter une fonction à tous les languages de programmation (unix/python/R etc.)
```

``` r
fnFs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"

``` r
## [1] "./MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"   "./MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
## [3] "./MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"
fnRs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"

``` r
## [1] "./MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"   "./MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
## [3] "./MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

##Infer sequence variants ###Dereplication

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz

    ## Encountered 1979 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz

    ## Encountered 1639 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz

    ## Encountered 1477 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz

    ## Encountered 904 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz

    ## Encountered 939 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz

    ## Encountered 1267 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz

    ## Encountered 1756 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz

    ## Encountered 1438 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz

    ## Encountered 3590 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz

    ## Encountered 2762 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz

    ## Encountered 3021 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz

    ## Encountered 1566 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz

    ## Encountered 3707 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz

    ## Encountered 1479 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz

    ## Encountered 1195 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz

    ## Encountered 1832 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz

    ## Encountered 1183 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz

    ## Encountered 1382 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz

    ## Encountered 1709 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz

    ## Encountered 897 unique sequences from 4314 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz

    ## Encountered 1660 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz

    ## Encountered 1349 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz

    ## Encountered 1335 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz

    ## Encountered 853 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz

    ## Encountered 880 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz

    ## Encountered 1286 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz

    ## Encountered 1803 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz

    ## Encountered 1265 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz

    ## Encountered 3414 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz

    ## Encountered 2522 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz

    ## Encountered 2771 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz

    ## Encountered 1415 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz

    ## Encountered 3290 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz

    ## Encountered 1390 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz

    ## Encountered 1134 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz

    ## Encountered 1635 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz

    ## Encountered 1084 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz

    ## Encountered 1161 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz

    ## Encountered 1502 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz

    ## Encountered 732 unique sequences from 4314 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
plotErrors(errF)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
plotErrors(errR)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

##Construct sequence table and remove chimeras

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
```

``` r
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  85 186   5   2

``` bash
cd ~
wget  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
```

    ## --2021-12-24 18:52:48--  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137283333 (131M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138.1_train_set.fa.gz.1’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 15.7M 8s
    ##     50K .......... .......... .......... .......... ..........  0% 9.79M 11s
    ##    100K .......... .......... .......... .......... ..........  0% 39.0M 8s
    ##    150K .......... .......... .......... .......... ..........  0% 13.5M 9s
    ##    200K .......... .......... .......... .......... ..........  0% 23.5M 8s
    ##    250K .......... .......... .......... .......... ..........  0% 12.9M 8s
    ##    300K .......... .......... .......... .......... ..........  0%  109M 7s
    ##    350K .......... .......... .......... .......... ..........  0% 16.6M 7s
    ##    400K .......... .......... .......... .......... ..........  0% 24.7M 7s
    ##    450K .......... .......... .......... .......... ..........  0% 19.9M 7s
    ##    500K .......... .......... .......... .......... ..........  0% 46.1M 7s
    ##    550K .......... .......... .......... .......... ..........  0% 16.2M 7s
    ##    600K .......... .......... .......... .......... ..........  0%  109M 6s
    ##    650K .......... .......... .......... .......... ..........  0% 17.0M 6s
    ##    700K .......... .......... .......... .......... ..........  0% 95.3M 6s
    ##    750K .......... .......... .......... .......... ..........  0% 13.5M 6s
    ##    800K .......... .......... .......... .......... ..........  0% 39.2M 6s
    ##    850K .......... .......... .......... .......... ..........  0% 24.1M 6s
    ##    900K .......... .......... .......... .......... ..........  0% 14.8M 6s
    ##    950K .......... .......... .......... .......... ..........  0% 62.8M 6s
    ##   1000K .......... .......... .......... .......... ..........  0% 19.0M 6s
    ##   1050K .......... .......... .......... .......... ..........  0% 84.4M 6s
    ##   1100K .......... .......... .......... .......... ..........  0% 14.7M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 70.7M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 15.6M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 81.5M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 18.3M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 33.8M 6s
    ##   1400K .......... .......... .......... .......... ..........  1% 12.3M 6s
    ##   1450K .......... .......... .......... .......... ..........  1%  106M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 11.5M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 96.0M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 13.5M 6s
    ##   1650K .......... .......... .......... .......... ..........  1%  113M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 16.4M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 83.9M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 18.3M 6s
    ##   1850K .......... .......... .......... .......... ..........  1% 63.6M 6s
    ##   1900K .......... .......... .......... .......... ..........  1% 6.13M 6s
    ##   1950K .......... .......... .......... .......... ..........  1% 43.3M 6s
    ##   2000K .......... .......... .......... .......... ..........  1%  107M 6s
    ##   2050K .......... .......... .......... .......... ..........  1% 29.5M 6s
    ##   2100K .......... .......... .......... .......... ..........  1% 28.8M 6s
    ##   2150K .......... .......... .......... .......... ..........  1% 14.1M 6s
    ##   2200K .......... .......... .......... .......... ..........  1% 13.3M 6s
    ##   2250K .......... .......... .......... .......... ..........  1% 34.5M 6s
    ##   2300K .......... .......... .......... .......... ..........  1% 18.9M 6s
    ##   2350K .......... .......... .......... .......... ..........  1% 13.5M 6s
    ##   2400K .......... .......... .......... .......... ..........  1% 13.1M 6s
    ##   2450K .......... .......... .......... .......... ..........  1% 92.4M 6s
    ##   2500K .......... .......... .......... .......... ..........  1% 15.4M 6s
    ##   2550K .......... .......... .......... .......... ..........  1% 78.0M 6s
    ##   2600K .......... .......... .......... .......... ..........  1% 15.0M 6s
    ##   2650K .......... .......... .......... .......... ..........  2% 13.4M 6s
    ##   2700K .......... .......... .......... .......... ..........  2% 16.1M 6s
    ##   2750K .......... .......... .......... .......... ..........  2% 22.9M 6s
    ##   2800K .......... .......... .......... .......... ..........  2% 29.3M 6s
    ##   2850K .......... .......... .......... .......... ..........  2% 14.5M 6s
    ##   2900K .......... .......... .......... .......... ..........  2% 13.7M 6s
    ##   2950K .......... .......... .......... .......... ..........  2% 22.3M 6s
    ##   3000K .......... .......... .......... .......... ..........  2% 26.9M 6s
    ##   3050K .......... .......... .......... .......... ..........  2% 13.1M 6s
    ##   3100K .......... .......... .......... .......... ..........  2% 13.6M 6s
    ##   3150K .......... .......... .......... .......... ..........  2% 12.9M 6s
    ##   3200K .......... .......... .......... .......... ..........  2% 76.1M 6s
    ##   3250K .......... .......... .......... .......... ..........  2% 14.2M 6s
    ##   3300K .......... .......... .......... .......... ..........  2% 41.4M 6s
    ##   3350K .......... .......... .......... .......... ..........  2% 12.0M 6s
    ##   3400K .......... .......... .......... .......... ..........  2% 34.3M 6s
    ##   3450K .......... .......... .......... .......... ..........  2% 14.0M 6s
    ##   3500K .......... .......... .......... .......... ..........  2% 8.59M 6s
    ##   3550K .......... .......... .......... .......... ..........  2% 10.5M 6s
    ##   3600K .......... .......... .......... .......... ..........  2% 11.7M 6s
    ##   3650K .......... .......... .......... .......... ..........  2% 98.4M 6s
    ##   3700K .......... .......... .......... .......... ..........  2% 13.5M 6s
    ##   3750K .......... .......... .......... .......... ..........  2% 13.0M 6s
    ##   3800K .......... .......... .......... .......... ..........  2% 14.0M 7s
    ##   3850K .......... .......... .......... .......... ..........  2% 78.5M 6s
    ##   3900K .......... .......... .......... .......... ..........  2% 8.69M 7s
    ##   3950K .......... .......... .......... .......... ..........  2% 11.1M 7s
    ##   4000K .......... .......... .......... .......... ..........  3% 14.0M 7s
    ##   4050K .......... .......... .......... .......... ..........  3% 13.4M 7s
    ##   4100K .......... .......... .......... .......... ..........  3% 91.7M 7s
    ##   4150K .......... .......... .......... .......... ..........  3% 12.6M 7s
    ##   4200K .......... .......... .......... .......... ..........  3% 13.6M 7s
    ##   4250K .......... .......... .......... .......... ..........  3% 17.9M 7s
    ##   4300K .......... .......... .......... .......... ..........  3% 14.8M 7s
    ##   4350K .......... .......... .......... .......... ..........  3% 16.5M 7s
    ##   4400K .......... .......... .......... .......... ..........  3% 76.6M 7s
    ##   4450K .......... .......... .......... .......... ..........  3% 14.8M 7s
    ##   4500K .......... .......... .......... .......... ..........  3% 99.6M 7s
    ##   4550K .......... .......... .......... .......... ..........  3% 14.3M 7s
    ##   4600K .......... .......... .......... .......... ..........  3%  114M 7s
    ##   4650K .......... .......... .......... .......... ..........  3% 17.4M 7s
    ##   4700K .......... .......... .......... .......... ..........  3% 73.0M 7s
    ##   4750K .......... .......... .......... .......... ..........  3% 13.9M 7s
    ##   4800K .......... .......... .......... .......... ..........  3% 31.1M 6s
    ##   4850K .......... .......... .......... .......... ..........  3% 20.0M 6s
    ##   4900K .......... .......... .......... .......... ..........  3% 17.7M 6s
    ##   4950K .......... .......... .......... .......... ..........  3% 63.3M 6s
    ##   5000K .......... .......... .......... .......... ..........  3% 17.0M 6s
    ##   5050K .......... .......... .......... .......... ..........  3% 32.5M 6s
    ##   5100K .......... .......... .......... .......... ..........  3% 13.4M 6s
    ##   5150K .......... .......... .......... .......... ..........  3% 91.5M 6s
    ##   5200K .......... .......... .......... .......... ..........  3% 18.0M 6s
    ##   5250K .......... .......... .......... .......... ..........  3% 85.0M 6s
    ##   5300K .......... .......... .......... .......... ..........  3% 16.0M 6s
    ##   5350K .......... .......... .......... .......... ..........  4% 67.5M 6s
    ##   5400K .......... .......... .......... .......... ..........  4% 15.0M 6s
    ##   5450K .......... .......... .......... .......... ..........  4% 78.5M 6s
    ##   5500K .......... .......... .......... .......... ..........  4% 18.4M 6s
    ##   5550K .......... .......... .......... .......... ..........  4% 56.9M 6s
    ##   5600K .......... .......... .......... .......... ..........  4% 11.2M 6s
    ##   5650K .......... .......... .......... .......... ..........  4% 99.8M 6s
    ##   5700K .......... .......... .......... .......... ..........  4% 12.2M 6s
    ##   5750K .......... .......... .......... .......... ..........  4% 85.5M 6s
    ##   5800K .......... .......... .......... .......... ..........  4% 12.3M 6s
    ##   5850K .......... .......... .......... .......... ..........  4%  103M 6s
    ##   5900K .......... .......... .......... .......... ..........  4% 15.4M 6s
    ##   5950K .......... .......... .......... .......... ..........  4% 66.5M 6s
    ##   6000K .......... .......... .......... .......... ..........  4% 13.0M 6s
    ##   6050K .......... .......... .......... .......... ..........  4% 90.6M 6s
    ##   6100K .......... .......... .......... .......... ..........  4% 14.1M 6s
    ##   6150K .......... .......... .......... .......... ..........  4% 78.5M 6s
    ##   6200K .......... .......... .......... .......... ..........  4% 13.8M 6s
    ##   6250K .......... .......... .......... .......... ..........  4%  102M 6s
    ##   6300K .......... .......... .......... .......... ..........  4% 11.7M 6s
    ##   6350K .......... .......... .......... .......... ..........  4% 91.3M 6s
    ##   6400K .......... .......... .......... .......... ..........  4% 12.6M 6s
    ##   6450K .......... .......... .......... .......... ..........  4% 11.1M 6s
    ##   6500K .......... .......... .......... .......... ..........  4% 85.4M 6s
    ##   6550K .......... .......... .......... .......... ..........  4% 93.7M 6s
    ##   6600K .......... .......... .......... .......... ..........  4% 13.0M 6s
    ##   6650K .......... .......... .......... .......... ..........  4% 13.8M 6s
    ##   6700K .......... .......... .......... .......... ..........  5%  146M 6s
    ##   6750K .......... .......... .......... .......... ..........  5% 12.1M 6s
    ##   6800K .......... .......... .......... .......... ..........  5%  157M 6s
    ##   6850K .......... .......... .......... .......... ..........  5% 13.0M 6s
    ##   6900K .......... .......... .......... .......... ..........  5% 88.8M 6s
    ##   6950K .......... .......... .......... .......... ..........  5% 19.1M 6s
    ##   7000K .......... .......... .......... .......... ..........  5% 30.1M 6s
    ##   7050K .......... .......... .......... .......... ..........  5% 22.1M 6s
    ##   7100K .......... .......... .......... .......... ..........  5% 28.4M 6s
    ##   7150K .......... .......... .......... .......... ..........  5% 23.5M 6s
    ##   7200K .......... .......... .......... .......... ..........  5% 40.2M 6s
    ##   7250K .......... .......... .......... .......... ..........  5%  117M 6s
    ##   7300K .......... .......... .......... .......... ..........  5% 13.1M 6s
    ##   7350K .......... .......... .......... .......... ..........  5% 42.5M 6s
    ##   7400K .......... .......... .......... .......... ..........  5% 20.8M 6s
    ##   7450K .......... .......... .......... .......... ..........  5%  113M 6s
    ##   7500K .......... .......... .......... .......... ..........  5% 15.4M 6s
    ##   7550K .......... .......... .......... .......... ..........  5%  107M 6s
    ##   7600K .......... .......... .......... .......... ..........  5% 15.1M 6s
    ##   7650K .......... .......... .......... .......... ..........  5%  135M 6s
    ##   7700K .......... .......... .......... .......... ..........  5% 15.2M 6s
    ##   7750K .......... .......... .......... .......... ..........  5%  154M 6s
    ##   7800K .......... .......... .......... .......... ..........  5% 14.9M 6s
    ##   7850K .......... .......... .......... .......... ..........  5% 75.2M 6s
    ##   7900K .......... .......... .......... .......... ..........  5% 16.8M 6s
    ##   7950K .......... .......... .......... .......... ..........  5%  104M 6s
    ##   8000K .......... .......... .......... .......... ..........  6% 13.7M 6s
    ##   8050K .......... .......... .......... .......... ..........  6%  121M 6s
    ##   8100K .......... .......... .......... .......... ..........  6% 14.1M 6s
    ##   8150K .......... .......... .......... .......... ..........  6%  112M 6s
    ##   8200K .......... .......... .......... .......... ..........  6% 11.6M 6s
    ##   8250K .......... .......... .......... .......... ..........  6%  164M 6s
    ##   8300K .......... .......... .......... .......... ..........  6% 11.6M 6s
    ##   8350K .......... .......... .......... .......... ..........  6%  143M 6s
    ##   8400K .......... .......... .......... .......... ..........  6% 11.3M 6s
    ##   8450K .......... .......... .......... .......... ..........  6%  155M 6s
    ##   8500K .......... .......... .......... .......... ..........  6% 13.0M 6s
    ##   8550K .......... .......... .......... .......... ..........  6%  107M 6s
    ##   8600K .......... .......... .......... .......... ..........  6% 13.4M 6s
    ##   8650K .......... .......... .......... .......... ..........  6%  121M 6s
    ##   8700K .......... .......... .......... .......... ..........  6%  152M 6s
    ##   8750K .......... .......... .......... .......... ..........  6% 12.2M 6s
    ##   8800K .......... .......... .......... .......... ..........  6% 11.2M 6s
    ##   8850K .......... .......... .......... .......... ..........  6%  102M 6s
    ##   8900K .......... .......... .......... .......... ..........  6%  103M 6s
    ##   8950K .......... .......... .......... .......... ..........  6% 17.1M 6s
    ##   9000K .......... .......... .......... .......... ..........  6% 55.0M 6s
    ##   9050K .......... .......... .......... .......... ..........  6%  126M 6s
    ##   9100K .......... .......... .......... .......... ..........  6% 18.4M 6s
    ##   9150K .......... .......... .......... .......... ..........  6% 50.2M 6s
    ##   9200K .......... .......... .......... .......... ..........  6% 16.7M 6s
    ##   9250K .......... .......... .......... .......... ..........  6% 65.5M 6s
    ##   9300K .......... .......... .......... .......... ..........  6%  104M 6s
    ##   9350K .......... .......... .......... .......... ..........  7% 12.7M 6s
    ##   9400K .......... .......... .......... .......... ..........  7%  163M 6s
    ##   9450K .......... .......... .......... .......... ..........  7% 12.5M 6s
    ##   9500K .......... .......... .......... .......... ..........  7%  128M 6s
    ##   9550K .......... .......... .......... .......... ..........  7%  134M 6s
    ##   9600K .......... .......... .......... .......... ..........  7% 16.3M 6s
    ##   9650K .......... .......... .......... .......... ..........  7% 72.5M 6s
    ##   9700K .......... .......... .......... .......... ..........  7% 24.2M 6s
    ##   9750K .......... .......... .......... .......... ..........  7% 28.3M 6s
    ##   9800K .......... .......... .......... .......... ..........  7%  104M 6s
    ##   9850K .......... .......... .......... .......... ..........  7% 15.2M 6s
    ##   9900K .......... .......... .......... .......... ..........  7%  128M 6s
    ##   9950K .......... .......... .......... .......... ..........  7% 14.9M 6s
    ##  10000K .......... .......... .......... .......... ..........  7%  102M 6s
    ##  10050K .......... .......... .......... .......... ..........  7%  172M 5s
    ##  10100K .......... .......... .......... .......... ..........  7% 14.7M 5s
    ##  10150K .......... .......... .......... .......... ..........  7%  125M 5s
    ##  10200K .......... .......... .......... .......... ..........  7%  166M 5s
    ##  10250K .......... .......... .......... .......... ..........  7% 15.1M 5s
    ##  10300K .......... .......... .......... .......... ..........  7%  147M 5s
    ##  10350K .......... .......... .......... .......... ..........  7% 18.4M 5s
    ##  10400K .......... .......... .......... .......... ..........  7% 50.6M 5s
    ##  10450K .......... .......... .......... .......... ..........  7%  132M 5s
    ##  10500K .......... .......... .......... .......... ..........  7% 19.4M 5s
    ##  10550K .......... .......... .......... .......... ..........  7% 49.3M 5s
    ##  10600K .......... .......... .......... .......... ..........  7% 21.1M 5s
    ##  10650K .......... .......... .......... .......... ..........  7% 39.2M 5s
    ##  10700K .......... .......... .......... .......... ..........  8%  173M 5s
    ##  10750K .......... .......... .......... .......... ..........  8% 20.9M 5s
    ##  10800K .......... .......... .......... .......... ..........  8% 42.4M 5s
    ##  10850K .......... .......... .......... .......... ..........  8% 21.3M 5s
    ##  10900K .......... .......... .......... .......... ..........  8% 44.6M 5s
    ##  10950K .......... .......... .......... .......... ..........  8%  136M 5s
    ##  11000K .......... .......... .......... .......... ..........  8% 23.7M 5s
    ##  11050K .......... .......... .......... .......... ..........  8% 37.4M 5s
    ##  11100K .......... .......... .......... .......... ..........  8% 34.6M 5s
    ##  11150K .......... .......... .......... .......... ..........  8% 27.0M 5s
    ##  11200K .......... .......... .......... .......... ..........  8% 79.4M 5s
    ##  11250K .......... .......... .......... .......... ..........  8% 29.0M 5s
    ##  11300K .......... .......... .......... .......... ..........  8% 39.9M 5s
    ##  11350K .......... .......... .......... .......... ..........  8% 29.1M 5s
    ##  11400K .......... .......... .......... .......... ..........  8% 54.9M 5s
    ##  11450K .......... .......... .......... .......... ..........  8% 20.6M 5s
    ##  11500K .......... .......... .......... .......... ..........  8% 33.5M 5s
    ##  11550K .......... .......... .......... .......... ..........  8% 26.4M 5s
    ##  11600K .......... .......... .......... .......... ..........  8% 25.8M 5s
    ##  11650K .......... .......... .......... .......... ..........  8%  148M 5s
    ##  11700K .......... .......... .......... .......... ..........  8% 27.4M 5s
    ##  11750K .......... .......... .......... .......... ..........  8% 24.2M 5s
    ##  11800K .......... .......... .......... .......... ..........  8% 26.5M 5s
    ##  11850K .......... .......... .......... .......... ..........  8%  168M 5s
    ##  11900K .......... .......... .......... .......... ..........  8% 32.6M 5s
    ##  11950K .......... .......... .......... .......... ..........  8% 27.8M 5s
    ##  12000K .......... .......... .......... .......... ..........  8% 22.7M 5s
    ##  12050K .......... .......... .......... .......... ..........  9% 28.9M 5s
    ##  12100K .......... .......... .......... .......... ..........  9%  147M 5s
    ##  12150K .......... .......... .......... .......... ..........  9% 28.3M 5s
    ##  12200K .......... .......... .......... .......... ..........  9% 31.5M 5s
    ##  12250K .......... .......... .......... .......... ..........  9% 38.8M 5s
    ##  12300K .......... .......... .......... .......... ..........  9%  140M 5s
    ##  12350K .......... .......... .......... .......... ..........  9% 23.7M 5s
    ##  12400K .......... .......... .......... .......... ..........  9% 34.7M 5s
    ##  12450K .......... .......... .......... .......... ..........  9% 24.9M 5s
    ##  12500K .......... .......... .......... .......... ..........  9% 40.4M 5s
    ##  12550K .......... .......... .......... .......... ..........  9%  143M 5s
    ##  12600K .......... .......... .......... .......... ..........  9% 24.8M 5s
    ##  12650K .......... .......... .......... .......... ..........  9% 35.4M 5s
    ##  12700K .......... .......... .......... .......... ..........  9% 43.8M 5s
    ##  12750K .......... .......... .......... .......... ..........  9% 12.4M 5s
    ##  12800K .......... .......... .......... .......... ..........  9%  132M 5s
    ##  12850K .......... .......... .......... .......... ..........  9%  178M 5s
    ##  12900K .......... .......... .......... .......... ..........  9% 11.2M 5s
    ##  12950K .......... .......... .......... .......... ..........  9%  145M 5s
    ##  13000K .......... .......... .......... .......... ..........  9% 16.6M 5s
    ##  13050K .......... .......... .......... .......... ..........  9% 73.4M 5s
    ##  13100K .......... .......... .......... .......... ..........  9%  133M 5s
    ##  13150K .......... .......... .......... .......... ..........  9% 17.0M 5s
    ##  13200K .......... .......... .......... .......... ..........  9% 95.8M 5s
    ##  13250K .......... .......... .......... .......... ..........  9% 17.8M 5s
    ##  13300K .......... .......... .......... .......... ..........  9%  133M 5s
    ##  13350K .......... .......... .......... .......... ..........  9%  113M 5s
    ##  13400K .......... .......... .......... .......... .......... 10% 18.2M 5s
    ##  13450K .......... .......... .......... .......... .......... 10% 70.5M 5s
    ##  13500K .......... .......... .......... .......... .......... 10%  142M 5s
    ##  13550K .......... .......... .......... .......... .......... 10% 19.5M 5s
    ##  13600K .......... .......... .......... .......... .......... 10% 44.3M 5s
    ##  13650K .......... .......... .......... .......... .......... 10% 38.9M 5s
    ##  13700K .......... .......... .......... .......... .......... 10% 28.9M 5s
    ##  13750K .......... .......... .......... .......... .......... 10%  127M 5s
    ##  13800K .......... .......... .......... .......... .......... 10% 35.2M 5s
    ##  13850K .......... .......... .......... .......... .......... 10% 24.8M 5s
    ##  13900K .......... .......... .......... .......... .......... 10% 34.3M 5s
    ##  13950K .......... .......... .......... .......... .......... 10% 23.2M 5s
    ##  14000K .......... .......... .......... .......... .......... 10%  146M 5s
    ##  14050K .......... .......... .......... .......... .......... 10% 36.5M 5s
    ##  14100K .......... .......... .......... .......... .......... 10% 25.1M 5s
    ##  14150K .......... .......... .......... .......... .......... 10% 49.9M 5s
    ##  14200K .......... .......... .......... .......... .......... 10%  127M 5s
    ##  14250K .......... .......... .......... .......... .......... 10% 23.9M 5s
    ##  14300K .......... .......... .......... .......... .......... 10% 39.1M 5s
    ##  14350K .......... .......... .......... .......... .......... 10% 21.6M 5s
    ##  14400K .......... .......... .......... .......... .......... 10% 60.9M 5s
    ##  14450K .......... .......... .......... .......... .......... 10%  127M 5s
    ##  14500K .......... .......... .......... .......... .......... 10% 23.0M 5s
    ##  14550K .......... .......... .......... .......... .......... 10% 50.5M 5s
    ##  14600K .......... .......... .......... .......... .......... 10% 32.9M 5s
    ##  14650K .......... .......... .......... .......... .......... 10% 66.6M 5s
    ##  14700K .......... .......... .......... .......... .......... 11% 35.1M 5s
    ##  14750K .......... .......... .......... .......... .......... 11% 10.8M 5s
    ##  14800K .......... .......... .......... .......... .......... 11%  147M 5s
    ##  14850K .......... .......... .......... .......... .......... 11% 54.8M 5s
    ##  14900K .......... .......... .......... .......... .......... 11% 20.4M 5s
    ##  14950K .......... .......... .......... .......... .......... 11%  134M 5s
    ##  15000K .......... .......... .......... .......... .......... 11% 32.7M 5s
    ##  15050K .......... .......... .......... .......... .......... 11% 22.1M 5s
    ##  15100K .......... .......... .......... .......... .......... 11% 57.7M 5s
    ##  15150K .......... .......... .......... .......... .......... 11% 18.9M 5s
    ##  15200K .......... .......... .......... .......... .......... 11%  149M 5s
    ##  15250K .......... .......... .......... .......... .......... 11% 36.8M 5s
    ##  15300K .......... .......... .......... .......... .......... 11% 21.8M 5s
    ##  15350K .......... .......... .......... .......... .......... 11% 32.5M 5s
    ##  15400K .......... .......... .......... .......... .......... 11%  166M 5s
    ##  15450K .......... .......... .......... .......... .......... 11% 22.2M 5s
    ##  15500K .......... .......... .......... .......... .......... 11% 45.6M 5s
    ##  15550K .......... .......... .......... .......... .......... 11% 20.8M 5s
    ##  15600K .......... .......... .......... .......... .......... 11% 33.7M 5s
    ##  15650K .......... .......... .......... .......... .......... 11%  149M 5s
    ##  15700K .......... .......... .......... .......... .......... 11% 22.0M 5s
    ##  15750K .......... .......... .......... .......... .......... 11% 35.9M 5s
    ##  15800K .......... .......... .......... .......... .......... 11% 21.8M 5s
    ##  15850K .......... .......... .......... .......... .......... 11%  154M 5s
    ##  15900K .......... .......... .......... .......... .......... 11% 35.0M 5s
    ##  15950K .......... .......... .......... .......... .......... 11% 28.0M 5s
    ##  16000K .......... .......... .......... .......... .......... 11% 26.5M 5s
    ##  16050K .......... .......... .......... .......... .......... 12% 49.8M 5s
    ##  16100K .......... .......... .......... .......... .......... 12%  101M 5s
    ##  16150K .......... .......... .......... .......... .......... 12% 23.2M 5s
    ##  16200K .......... .......... .......... .......... .......... 12% 27.7M 5s
    ##  16250K .......... .......... .......... .......... .......... 12% 25.0M 5s
    ##  16300K .......... .......... .......... .......... .......... 12%  160M 5s
    ##  16350K .......... .......... .......... .......... .......... 12% 28.3M 5s
    ##  16400K .......... .......... .......... .......... .......... 12% 34.0M 5s
    ##  16450K .......... .......... .......... .......... .......... 12% 15.7M 5s
    ##  16500K .......... .......... .......... .......... .......... 12% 61.1M 4s
    ##  16550K .......... .......... .......... .......... .......... 12% 85.9M 4s
    ##  16600K .......... .......... .......... .......... .......... 12% 18.7M 4s
    ##  16650K .......... .......... .......... .......... .......... 12% 71.1M 4s
    ##  16700K .......... .......... .......... .......... .......... 12% 19.3M 4s
    ##  16750K .......... .......... .......... .......... .......... 12% 70.9M 4s
    ##  16800K .......... .......... .......... .......... .......... 12% 96.4M 4s
    ##  16850K .......... .......... .......... .......... .......... 12% 19.3M 4s
    ##  16900K .......... .......... .......... .......... .......... 12% 76.0M 4s
    ##  16950K .......... .......... .......... .......... .......... 12% 20.0M 4s
    ##  17000K .......... .......... .......... .......... .......... 12% 45.2M 4s
    ##  17050K .......... .......... .......... .......... .......... 12% 86.4M 4s
    ##  17100K .......... .......... .......... .......... .......... 12% 22.0M 4s
    ##  17150K .......... .......... .......... .......... .......... 12% 68.5M 4s
    ##  17200K .......... .......... .......... .......... .......... 12% 18.3M 4s
    ##  17250K .......... .......... .......... .......... .......... 12% 84.3M 4s
    ##  17300K .......... .......... .......... .......... .......... 12% 96.5M 4s
    ##  17350K .......... .......... .......... .......... .......... 12% 17.1M 4s
    ##  17400K .......... .......... .......... .......... .......... 13% 66.1M 4s
    ##  17450K .......... .......... .......... .......... .......... 13% 14.9M 4s
    ##  17500K .......... .......... .......... .......... .......... 13% 90.1M 4s
    ##  17550K .......... .......... .......... .......... .......... 13% 89.0M 4s
    ##  17600K .......... .......... .......... .......... .......... 13% 17.2M 4s
    ##  17650K .......... .......... .......... .......... .......... 13% 98.4M 4s
    ##  17700K .......... .......... .......... .......... .......... 13% 21.9M 4s
    ##  17750K .......... .......... .......... .......... .......... 13% 61.4M 4s
    ##  17800K .......... .......... .......... .......... .......... 13%  111M 4s
    ##  17850K .......... .......... .......... .......... .......... 13% 16.7M 4s
    ##  17900K .......... .......... .......... .......... .......... 13% 90.5M 4s
    ##  17950K .......... .......... .......... .......... .......... 13% 15.4M 4s
    ##  18000K .......... .......... .......... .......... .......... 13% 80.0M 4s
    ##  18050K .......... .......... .......... .......... .......... 13%  101M 4s
    ##  18100K .......... .......... .......... .......... .......... 13% 99.7M 4s
    ##  18150K .......... .......... .......... .......... .......... 13% 18.6M 4s
    ##  18200K .......... .......... .......... .......... .......... 13% 94.6M 4s
    ##  18250K .......... .......... .......... .......... .......... 13% 60.0M 4s
    ##  18300K .......... .......... .......... .......... .......... 13% 45.5M 4s
    ##  18350K .......... .......... .......... .......... .......... 13% 22.0M 4s
    ##  18400K .......... .......... .......... .......... .......... 13% 36.3M 4s
    ##  18450K .......... .......... .......... .......... .......... 13% 97.2M 4s
    ##  18500K .......... .......... .......... .......... .......... 13% 23.3M 4s
    ##  18550K .......... .......... .......... .......... .......... 13% 26.9M 4s
    ##  18600K .......... .......... .......... .......... .......... 13% 37.3M 4s
    ##  18650K .......... .......... .......... .......... .......... 13% 77.2M 4s
    ##  18700K .......... .......... .......... .......... .......... 13% 20.4M 4s
    ##  18750K .......... .......... .......... .......... .......... 14% 39.0M 4s
    ##  18800K .......... .......... .......... .......... .......... 14% 20.0M 4s
    ##  18850K .......... .......... .......... .......... .......... 14% 79.1M 4s
    ##  18900K .......... .......... .......... .......... .......... 14% 92.9M 4s
    ##  18950K .......... .......... .......... .......... .......... 14% 18.7M 4s
    ##  19000K .......... .......... .......... .......... .......... 14% 91.3M 4s
    ##  19050K .......... .......... .......... .......... .......... 14% 18.3M 4s
    ##  19100K .......... .......... .......... .......... .......... 14% 74.3M 4s
    ##  19150K .......... .......... .......... .......... .......... 14% 15.9M 4s
    ##  19200K .......... .......... .......... .......... .......... 14% 55.2M 4s
    ##  19250K .......... .......... .......... .......... .......... 14% 89.6M 4s
    ##  19300K .......... .......... .......... .......... .......... 14% 14.4M 4s
    ##  19350K .......... .......... .......... .......... .......... 14% 75.1M 4s
    ##  19400K .......... .......... .......... .......... .......... 14%  110M 4s
    ##  19450K .......... .......... .......... .......... .......... 14% 14.5M 4s
    ##  19500K .......... .......... .......... .......... .......... 14% 83.3M 4s
    ##  19550K .......... .......... .......... .......... .......... 14% 18.0M 4s
    ##  19600K .......... .......... .......... .......... .......... 14% 51.7M 4s
    ##  19650K .......... .......... .......... .......... .......... 14% 70.9M 4s
    ##  19700K .......... .......... .......... .......... .......... 14% 20.5M 4s
    ##  19750K .......... .......... .......... .......... .......... 14% 41.2M 4s
    ##  19800K .......... .......... .......... .......... .......... 14% 18.5M 4s
    ##  19850K .......... .......... .......... .......... .......... 14%  103M 4s
    ##  19900K .......... .......... .......... .......... .......... 14% 57.6M 4s
    ##  19950K .......... .......... .......... .......... .......... 14% 18.4M 4s
    ##  20000K .......... .......... .......... .......... .......... 14% 38.5M 4s
    ##  20050K .......... .......... .......... .......... .......... 14% 20.3M 4s
    ##  20100K .......... .......... .......... .......... .......... 15% 93.8M 4s
    ##  20150K .......... .......... .......... .......... .......... 15% 43.5M 4s
    ##  20200K .......... .......... .......... .......... .......... 15% 15.9M 4s
    ##  20250K .......... .......... .......... .......... .......... 15% 76.2M 4s
    ##  20300K .......... .......... .......... .......... .......... 15%  117M 4s
    ##  20350K .......... .......... .......... .......... .......... 15% 26.8M 4s
    ##  20400K .......... .......... .......... .......... .......... 15% 25.1M 4s
    ##  20450K .......... .......... .......... .......... .......... 15% 42.4M 4s
    ##  20500K .......... .......... .......... .......... .......... 15% 18.5M 4s
    ##  20550K .......... .......... .......... .......... .......... 15% 88.4M 4s
    ##  20600K .......... .......... .......... .......... .......... 15% 51.0M 4s
    ##  20650K .......... .......... .......... .......... .......... 15% 14.5M 4s
    ##  20700K .......... .......... .......... .......... .......... 15% 80.6M 4s
    ##  20750K .......... .......... .......... .......... .......... 15% 19.5M 4s
    ##  20800K .......... .......... .......... .......... .......... 15% 72.4M 4s
    ##  20850K .......... .......... .......... .......... .......... 15% 89.4M 4s
    ##  20900K .......... .......... .......... .......... .......... 15% 18.9M 4s
    ##  20950K .......... .......... .......... .......... .......... 15% 73.8M 4s
    ##  21000K .......... .......... .......... .......... .......... 15%  100M 4s
    ##  21050K .......... .......... .......... .......... .......... 15% 20.2M 4s
    ##  21100K .......... .......... .......... .......... .......... 15% 81.4M 4s
    ##  21150K .......... .......... .......... .......... .......... 15% 20.9M 4s
    ##  21200K .......... .......... .......... .......... .......... 15% 40.2M 4s
    ##  21250K .......... .......... .......... .......... .......... 15% 24.2M 4s
    ##  21300K .......... .......... .......... .......... .......... 15% 95.3M 4s
    ##  21350K .......... .......... .......... .......... .......... 15% 22.4M 4s
    ##  21400K .......... .......... .......... .......... .......... 15% 96.9M 4s
    ##  21450K .......... .......... .......... .......... .......... 16% 14.1M 4s
    ##  21500K .......... .......... .......... .......... .......... 16% 71.3M 4s
    ##  21550K .......... .......... .......... .......... .......... 16% 91.8M 4s
    ##  21600K .......... .......... .......... .......... .......... 16% 14.4M 4s
    ##  21650K .......... .......... .......... .......... .......... 16% 92.2M 4s
    ##  21700K .......... .......... .......... .......... .......... 16% 14.1M 4s
    ##  21750K .......... .......... .......... .......... .......... 16% 75.0M 4s
    ##  21800K .......... .......... .......... .......... .......... 16%  103M 4s
    ##  21850K .......... .......... .......... .......... .......... 16% 16.2M 4s
    ##  21900K .......... .......... .......... .......... .......... 16% 82.2M 4s
    ##  21950K .......... .......... .......... .......... .......... 16% 17.6M 4s
    ##  22000K .......... .......... .......... .......... .......... 16% 54.9M 4s
    ##  22050K .......... .......... .......... .......... .......... 16% 73.0M 4s
    ##  22100K .......... .......... .......... .......... .......... 16% 20.2M 4s
    ##  22150K .......... .......... .......... .......... .......... 16% 59.4M 4s
    ##  22200K .......... .......... .......... .......... .......... 16% 77.6M 4s
    ##  22250K .......... .......... .......... .......... .......... 16% 17.0M 4s
    ##  22300K .......... .......... .......... .......... .......... 16% 72.6M 4s
    ##  22350K .......... .......... .......... .......... .......... 16% 15.3M 4s
    ##  22400K .......... .......... .......... .......... .......... 16% 69.2M 4s
    ##  22450K .......... .......... .......... .......... .......... 16%  113M 4s
    ##  22500K .......... .......... .......... .......... .......... 16% 13.6M 4s
    ##  22550K .......... .......... .......... .......... .......... 16% 99.8M 4s
    ##  22600K .......... .......... .......... .......... .......... 16% 15.1M 4s
    ##  22650K .......... .......... .......... .......... .......... 16% 77.8M 4s
    ##  22700K .......... .......... .......... .......... .......... 16% 63.6M 4s
    ##  22750K .......... .......... .......... .......... .......... 17% 72.4M 4s
    ##  22800K .......... .......... .......... .......... .......... 17% 22.9M 4s
    ##  22850K .......... .......... .......... .......... .......... 17% 46.8M 4s
    ##  22900K .......... .......... .......... .......... .......... 17%  106M 4s
    ##  22950K .......... .......... .......... .......... .......... 17% 20.9M 4s
    ##  23000K .......... .......... .......... .......... .......... 17% 39.9M 4s
    ##  23050K .......... .......... .......... .......... .......... 17% 75.3M 4s
    ##  23100K .......... .......... .......... .......... .......... 17%  116M 4s
    ##  23150K .......... .......... .......... .......... .......... 17% 38.0M 4s
    ##  23200K .......... .......... .......... .......... .......... 17% 22.4M 4s
    ##  23250K .......... .......... .......... .......... .......... 17% 84.3M 4s
    ##  23300K .......... .......... .......... .......... .......... 17% 42.4M 4s
    ##  23350K .......... .......... .......... .......... .......... 17% 30.9M 4s
    ##  23400K .......... .......... .......... .......... .......... 17% 85.0M 4s
    ##  23450K .......... .......... .......... .......... .......... 17% 25.6M 4s
    ##  23500K .......... .......... .......... .......... .......... 17% 82.9M 4s
    ##  23550K .......... .......... .......... .......... .......... 17% 62.8M 4s
    ##  23600K .......... .......... .......... .......... .......... 17% 22.7M 4s
    ##  23650K .......... .......... .......... .......... .......... 17% 93.3M 4s
    ##  23700K .......... .......... .......... .......... .......... 17% 68.6M 4s
    ##  23750K .......... .......... .......... .......... .......... 17% 39.1M 4s
    ##  23800K .......... .......... .......... .......... .......... 17% 30.1M 4s
    ##  23850K .......... .......... .......... .......... .......... 17% 66.2M 4s
    ##  23900K .......... .......... .......... .......... .......... 17% 74.7M 4s
    ##  23950K .......... .......... .......... .......... .......... 17% 44.9M 4s
    ##  24000K .......... .......... .......... .......... .......... 17% 21.2M 4s
    ##  24050K .......... .......... .......... .......... .......... 17% 39.4M 4s
    ##  24100K .......... .......... .......... .......... .......... 18% 81.8M 4s
    ##  24150K .......... .......... .......... .......... .......... 18%  105M 4s
    ##  24200K .......... .......... .......... .......... .......... 18% 41.2M 4s
    ##  24250K .......... .......... .......... .......... .......... 18% 49.6M 4s
    ##  24300K .......... .......... .......... .......... .......... 18% 55.1M 4s
    ##  24350K .......... .......... .......... .......... .......... 18% 51.9M 4s
    ##  24400K .......... .......... .......... .......... .......... 18% 55.8M 4s
    ##  24450K .......... .......... .......... .......... .......... 18% 28.7M 4s
    ##  24500K .......... .......... .......... .......... .......... 18% 42.3M 4s
    ##  24550K .......... .......... .......... .......... .......... 18% 82.8M 4s
    ##  24600K .......... .......... .......... .......... .......... 18% 67.6M 4s
    ##  24650K .......... .......... .......... .......... .......... 18% 22.6M 4s
    ##  24700K .......... .......... .......... .......... .......... 18% 61.3M 4s
    ##  24750K .......... .......... .......... .......... .......... 18% 41.4M 4s
    ##  24800K .......... .......... .......... .......... .......... 18% 29.1M 4s
    ##  24850K .......... .......... .......... .......... .......... 18% 96.7M 4s
    ##  24900K .......... .......... .......... .......... .......... 18% 53.2M 4s
    ##  24950K .......... .......... .......... .......... .......... 18% 29.3M 4s
    ##  25000K .......... .......... .......... .......... .......... 18% 46.1M 4s
    ##  25050K .......... .......... .......... .......... .......... 18% 63.7M 4s
    ##  25100K .......... .......... .......... .......... .......... 18%  111M 4s
    ##  25150K .......... .......... .......... .......... .......... 18% 11.3M 4s
    ##  25200K .......... .......... .......... .......... .......... 18% 97.6M 4s
    ##  25250K .......... .......... .......... .......... .......... 18%  113M 4s
    ##  25300K .......... .......... .......... .......... .......... 18%  101M 4s
    ##  25350K .......... .......... .......... .......... .......... 18% 17.9M 4s
    ##  25400K .......... .......... .......... .......... .......... 18% 51.1M 4s
    ##  25450K .......... .......... .......... .......... .......... 19% 59.1M 4s
    ##  25500K .......... .......... .......... .......... .......... 19% 89.7M 4s
    ##  25550K .......... .......... .......... .......... .......... 19% 27.9M 4s
    ##  25600K .......... .......... .......... .......... .......... 19% 30.5M 4s
    ##  25650K .......... .......... .......... .......... .......... 19%  105M 4s
    ##  25700K .......... .......... .......... .......... .......... 19%  117M 4s
    ##  25750K .......... .......... .......... .......... .......... 19% 43.0M 4s
    ##  25800K .......... .......... .......... .......... .......... 19% 14.9M 4s
    ##  25850K .......... .......... .......... .......... .......... 19% 91.5M 4s
    ##  25900K .......... .......... .......... .......... .......... 19% 84.3M 4s
    ##  25950K .......... .......... .......... .......... .......... 19% 19.4M 4s
    ##  26000K .......... .......... .......... .......... .......... 19% 71.2M 4s
    ##  26050K .......... .......... .......... .......... .......... 19% 79.4M 4s
    ##  26100K .......... .......... .......... .......... .......... 19% 71.7M 4s
    ##  26150K .......... .......... .......... .......... .......... 19% 76.0M 4s
    ##  26200K .......... .......... .......... .......... .......... 19% 25.4M 4s
    ##  26250K .......... .......... .......... .......... .......... 19% 36.2M 4s
    ##  26300K .......... .......... .......... .......... .......... 19% 87.7M 4s
    ##  26350K .......... .......... .......... .......... .......... 19% 68.6M 4s
    ##  26400K .......... .......... .......... .......... .......... 19% 18.6M 4s
    ##  26450K .......... .......... .......... .......... .......... 19% 82.9M 4s
    ##  26500K .......... .......... .......... .......... .......... 19% 94.9M 4s
    ##  26550K .......... .......... .......... .......... .......... 19% 70.4M 4s
    ##  26600K .......... .......... .......... .......... .......... 19% 19.5M 4s
    ##  26650K .......... .......... .......... .......... .......... 19% 77.2M 4s
    ##  26700K .......... .......... .......... .......... .......... 19%  104M 4s
    ##  26750K .......... .......... .......... .......... .......... 19% 54.6M 4s
    ##  26800K .......... .......... .......... .......... .......... 20% 20.9M 4s
    ##  26850K .......... .......... .......... .......... .......... 20% 81.3M 4s
    ##  26900K .......... .......... .......... .......... .......... 20%  117M 4s
    ##  26950K .......... .......... .......... .......... .......... 20% 95.1M 4s
    ##  27000K .......... .......... .......... .......... .......... 20% 23.1M 4s
    ##  27050K .......... .......... .......... .......... .......... 20% 90.0M 4s
    ##  27100K .......... .......... .......... .......... .......... 20% 98.6M 4s
    ##  27150K .......... .......... .......... .......... .......... 20% 90.8M 4s
    ##  27200K .......... .......... .......... .......... .......... 20% 16.8M 4s
    ##  27250K .......... .......... .......... .......... .......... 20%  102M 4s
    ##  27300K .......... .......... .......... .......... .......... 20%  104M 4s
    ##  27350K .......... .......... .......... .......... .......... 20%  100M 4s
    ##  27400K .......... .......... .......... .......... .......... 20% 22.4M 4s
    ##  27450K .......... .......... .......... .......... .......... 20% 91.5M 4s
    ##  27500K .......... .......... .......... .......... .......... 20% 86.9M 4s
    ##  27550K .......... .......... .......... .......... .......... 20% 19.0M 4s
    ##  27600K .......... .......... .......... .......... .......... 20% 69.5M 4s
    ##  27650K .......... .......... .......... .......... .......... 20% 92.3M 4s
    ##  27700K .......... .......... .......... .......... .......... 20% 86.9M 4s
    ##  27750K .......... .......... .......... .......... .......... 20% 21.9M 4s
    ##  27800K .......... .......... .......... .......... .......... 20% 69.2M 4s
    ##  27850K .......... .......... .......... .......... .......... 20% 92.0M 4s
    ##  27900K .......... .......... .......... .......... .......... 20% 65.3M 4s
    ##  27950K .......... .......... .......... .......... .......... 20% 24.6M 4s
    ##  28000K .......... .......... .......... .......... .......... 20% 47.6M 4s
    ##  28050K .......... .......... .......... .......... .......... 20% 80.7M 4s
    ##  28100K .......... .......... .......... .......... .......... 20%  122M 4s
    ##  28150K .......... .......... .......... .......... .......... 21% 22.9M 4s
    ##  28200K .......... .......... .......... .......... .......... 21% 72.4M 4s
    ##  28250K .......... .......... .......... .......... .......... 21% 72.1M 4s
    ##  28300K .......... .......... .......... .......... .......... 21%  123M 4s
    ##  28350K .......... .......... .......... .......... .......... 21% 13.5M 4s
    ##  28400K .......... .......... .......... .......... .......... 21%  120M 4s
    ##  28450K .......... .......... .......... .......... .......... 21%  129M 4s
    ##  28500K .......... .......... .......... .......... .......... 21% 13.6M 4s
    ##  28550K .......... .......... .......... .......... .......... 21% 67.0M 4s
    ##  28600K .......... .......... .......... .......... .......... 21% 68.5M 3s
    ##  28650K .......... .......... .......... .......... .......... 21%  122M 3s
    ##  28700K .......... .......... .......... .......... .......... 21% 23.9M 3s
    ##  28750K .......... .......... .......... .......... .......... 21% 67.8M 3s
    ##  28800K .......... .......... .......... .......... .......... 21%  107M 3s
    ##  28850K .......... .......... .......... .......... .......... 21%  102M 3s
    ##  28900K .......... .......... .......... .......... .......... 21% 26.2M 3s
    ##  28950K .......... .......... .......... .......... .......... 21% 81.6M 3s
    ##  29000K .......... .......... .......... .......... .......... 21% 98.2M 3s
    ##  29050K .......... .......... .......... .......... .......... 21% 67.2M 3s
    ##  29100K .......... .......... .......... .......... .......... 21% 33.0M 3s
    ##  29150K .......... .......... .......... .......... .......... 21% 25.6M 3s
    ##  29200K .......... .......... .......... .......... .......... 21% 83.2M 3s
    ##  29250K .......... .......... .......... .......... .......... 21%  117M 3s
    ##  29300K .......... .......... .......... .......... .......... 21% 40.2M 3s
    ##  29350K .......... .......... .......... .......... .......... 21% 47.4M 3s
    ##  29400K .......... .......... .......... .......... .......... 21% 43.6M 3s
    ##  29450K .......... .......... .......... .......... .......... 22%  121M 3s
    ##  29500K .......... .......... .......... .......... .......... 22% 54.2M 3s
    ##  29550K .......... .......... .......... .......... .......... 22% 37.4M 3s
    ##  29600K .......... .......... .......... .......... .......... 22% 67.5M 3s
    ##  29650K .......... .......... .......... .......... .......... 22% 64.7M 3s
    ##  29700K .......... .......... .......... .......... .......... 22% 81.9M 3s
    ##  29750K .......... .......... .......... .......... .......... 22% 30.8M 3s
    ##  29800K .......... .......... .......... .......... .......... 22% 40.5M 3s
    ##  29850K .......... .......... .......... .......... .......... 22% 95.9M 3s
    ##  29900K .......... .......... .......... .......... .......... 22% 28.0M 3s
    ##  29950K .......... .......... .......... .......... .......... 22% 27.4M 3s
    ##  30000K .......... .......... .......... .......... .......... 22% 86.0M 3s
    ##  30050K .......... .......... .......... .......... .......... 22% 59.1M 3s
    ##  30100K .......... .......... .......... .......... .......... 22%  100M 3s
    ##  30150K .......... .......... .......... .......... .......... 22% 20.0M 3s
    ##  30200K .......... .......... .......... .......... .......... 22% 76.8M 3s
    ##  30250K .......... .......... .......... .......... .......... 22%  111M 3s
    ##  30300K .......... .......... .......... .......... .......... 22%  130M 3s
    ##  30350K .......... .......... .......... .......... .......... 22% 21.5M 3s
    ##  30400K .......... .......... .......... .......... .......... 22%  107M 3s
    ##  30450K .......... .......... .......... .......... .......... 22%  101M 3s
    ##  30500K .......... .......... .......... .......... .......... 22%  130M 3s
    ##  30550K .......... .......... .......... .......... .......... 22% 21.7M 3s
    ##  30600K .......... .......... .......... .......... .......... 22% 90.7M 3s
    ##  30650K .......... .......... .......... .......... .......... 22%  114M 3s
    ##  30700K .......... .......... .......... .......... .......... 22%  113M 3s
    ##  30750K .......... .......... .......... .......... .......... 22% 23.2M 3s
    ##  30800K .......... .......... .......... .......... .......... 23% 93.1M 3s
    ##  30850K .......... .......... .......... .......... .......... 23%  103M 3s
    ##  30900K .......... .......... .......... .......... .......... 23%  127M 3s
    ##  30950K .......... .......... .......... .......... .......... 23% 23.5M 3s
    ##  31000K .......... .......... .......... .......... .......... 23% 86.7M 3s
    ##  31050K .......... .......... .......... .......... .......... 23% 86.9M 3s
    ##  31100K .......... .......... .......... .......... .......... 23%  122M 3s
    ##  31150K .......... .......... .......... .......... .......... 23% 31.9M 3s
    ##  31200K .......... .......... .......... .......... .......... 23% 24.3M 3s
    ##  31250K .......... .......... .......... .......... .......... 23%  106M 3s
    ##  31300K .......... .......... .......... .......... .......... 23% 45.0M 3s
    ##  31350K .......... .......... .......... .......... .......... 23%  106M 3s
    ##  31400K .......... .......... .......... .......... .......... 23% 34.4M 3s
    ##  31450K .......... .......... .......... .......... .......... 23%  124M 3s
    ##  31500K .......... .......... .......... .......... .......... 23% 36.1M 3s
    ##  31550K .......... .......... .......... .......... .......... 23% 55.7M 3s
    ##  31600K .......... .......... .......... .......... .......... 23% 71.2M 3s
    ##  31650K .......... .......... .......... .......... .......... 23%  100M 3s
    ##  31700K .......... .......... .......... .......... .......... 23% 23.3M 3s
    ##  31750K .......... .......... .......... .......... .......... 23% 30.0M 3s
    ##  31800K .......... .......... .......... .......... .......... 23% 89.5M 3s
    ##  31850K .......... .......... .......... .......... .......... 23%  109M 3s
    ##  31900K .......... .......... .......... .......... .......... 23% 27.2M 3s
    ##  31950K .......... .......... .......... .......... .......... 23% 28.9M 3s
    ##  32000K .......... .......... .......... .......... .......... 23%  114M 3s
    ##  32050K .......... .......... .......... .......... .......... 23% 31.9M 3s
    ##  32100K .......... .......... .......... .......... .......... 23% 93.5M 3s
    ##  32150K .......... .......... .......... .......... .......... 24% 20.9M 3s
    ##  32200K .......... .......... .......... .......... .......... 24% 97.4M 3s
    ##  32250K .......... .......... .......... .......... .......... 24% 97.1M 3s
    ##  32300K .......... .......... .......... .......... .......... 24% 6.28M 3s
    ##  32350K .......... .......... .......... .......... .......... 24% 84.3M 3s
    ##  32400K .......... .......... .......... .......... .......... 24% 95.3M 3s
    ##  32450K .......... .......... .......... .......... .......... 24% 89.2M 3s
    ##  32500K .......... .......... .......... .......... .......... 24% 91.3M 3s
    ##  32550K .......... .......... .......... .......... .......... 24% 92.2M 3s
    ##  32600K .......... .......... .......... .......... .......... 24% 96.5M 3s
    ##  32650K .......... .......... .......... .......... .......... 24%  107M 3s
    ##  32700K .......... .......... .......... .......... .......... 24%  126M 3s
    ##  32750K .......... .......... .......... .......... .......... 24% 26.6M 3s
    ##  32800K .......... .......... .......... .......... .......... 24% 99.8M 3s
    ##  32850K .......... .......... .......... .......... .......... 24%  123M 3s
    ##  32900K .......... .......... .......... .......... .......... 24% 14.3M 3s
    ##  32950K .......... .......... .......... .......... .......... 24% 75.3M 3s
    ##  33000K .......... .......... .......... .......... .......... 24%  127M 3s
    ##  33050K .......... .......... .......... .......... .......... 24% 18.0M 3s
    ##  33100K .......... .......... .......... .......... .......... 24%  115M 3s
    ##  33150K .......... .......... .......... .......... .......... 24% 17.1M 3s
    ##  33200K .......... .......... .......... .......... .......... 24%  108M 3s
    ##  33250K .......... .......... .......... .......... .......... 24%  128M 3s
    ##  33300K .......... .......... .......... .......... .......... 24% 16.7M 3s
    ##  33350K .......... .......... .......... .......... .......... 24% 52.6M 3s
    ##  33400K .......... .......... .......... .......... .......... 24% 17.1M 3s
    ##  33450K .......... .......... .......... .......... .......... 24%  119M 3s
    ##  33500K .......... .......... .......... .......... .......... 25% 81.7M 3s
    ##  33550K .......... .......... .......... .......... .......... 25% 10.5M 3s
    ##  33600K .......... .......... .......... .......... .......... 25%  124M 3s
    ##  33650K .......... .......... .......... .......... .......... 25% 13.5M 3s
    ##  33700K .......... .......... .......... .......... .......... 25% 99.2M 3s
    ##  33750K .......... .......... .......... .......... .......... 25% 86.0M 3s
    ##  33800K .......... .......... .......... .......... .......... 25% 51.0M 3s
    ##  33850K .......... .......... .......... .......... .......... 25% 17.1M 3s
    ##  33900K .......... .......... .......... .......... .......... 25%  125M 3s
    ##  33950K .......... .......... .......... .......... .......... 25% 15.4M 3s
    ##  34000K .......... .......... .......... .......... .......... 25% 80.3M 3s
    ##  34050K .......... .......... .......... .......... .......... 25% 5.80M 3s
    ##  34100K .......... .......... .......... .......... .......... 25%  113M 3s
    ##  34150K .......... .......... .......... .......... .......... 25%  120M 3s
    ##  34200K .......... .......... .......... .......... .......... 25%  139M 3s
    ##  34250K .......... .......... .......... .......... .......... 25%  122M 3s
    ##  34300K .......... .......... .......... .......... .......... 25%  138M 3s
    ##  34350K .......... .......... .......... .......... .......... 25% 26.3M 3s
    ##  34400K .......... .......... .......... .......... .......... 25% 24.3M 3s
    ##  34450K .......... .......... .......... .......... .......... 25% 92.0M 3s
    ##  34500K .......... .......... .......... .......... .......... 25% 43.4M 3s
    ##  34550K .......... .......... .......... .......... .......... 25% 17.0M 3s
    ##  34600K .......... .......... .......... .......... .......... 25% 14.3M 3s
    ##  34650K .......... .......... .......... .......... .......... 25% 93.8M 3s
    ##  34700K .......... .......... .......... .......... .......... 25%  130M 3s
    ##  34750K .......... .......... .......... .......... .......... 25% 13.4M 3s
    ##  34800K .......... .......... .......... .......... .......... 25% 11.6M 3s
    ##  34850K .......... .......... .......... .......... .......... 26%  128M 3s
    ##  34900K .......... .......... .......... .......... .......... 26% 14.5M 3s
    ##  34950K .......... .......... .......... .......... .......... 26% 71.1M 3s
    ##  35000K .......... .......... .......... .......... .......... 26%  121M 3s
    ##  35050K .......... .......... .......... .......... .......... 26% 14.9M 3s
    ##  35100K .......... .......... .......... .......... .......... 26%  110M 3s
    ##  35150K .......... .......... .......... .......... .......... 26% 14.5M 3s
    ##  35200K .......... .......... .......... .......... .......... 26% 98.7M 3s
    ##  35250K .......... .......... .......... .......... .......... 26% 14.5M 3s
    ##  35300K .......... .......... .......... .......... .......... 26%  108M 3s
    ##  35350K .......... .......... .......... .......... .......... 26% 19.5M 3s
    ##  35400K .......... .......... .......... .......... .......... 26% 31.5M 3s
    ##  35450K .......... .......... .......... .......... .......... 26% 35.3M 3s
    ##  35500K .......... .......... .......... .......... .......... 26% 25.1M 3s
    ##  35550K .......... .......... .......... .......... .......... 26% 18.0M 3s
    ##  35600K .......... .......... .......... .......... .......... 26% 27.9M 3s
    ##  35650K .......... .......... .......... .......... .......... 26% 87.5M 3s
    ##  35700K .......... .......... .......... .......... .......... 26% 13.1M 3s
    ##  35750K .......... .......... .......... .......... .......... 26% 30.1M 3s
    ##  35800K .......... .......... .......... .......... .......... 26% 19.2M 3s
    ##  35850K .......... .......... .......... .......... .......... 26% 87.1M 3s
    ##  35900K .......... .......... .......... .......... .......... 26% 30.4M 3s
    ##  35950K .......... .......... .......... .......... .......... 26% 25.7M 3s
    ##  36000K .......... .......... .......... .......... .......... 26% 19.4M 3s
    ##  36050K .......... .......... .......... .......... .......... 26% 83.8M 3s
    ##  36100K .......... .......... .......... .......... .......... 26% 17.9M 3s
    ##  36150K .......... .......... .......... .......... .......... 27% 36.2M 3s
    ##  36200K .......... .......... .......... .......... .......... 27% 13.0M 3s
    ##  36250K .......... .......... .......... .......... .......... 27% 98.2M 3s
    ##  36300K .......... .......... .......... .......... .......... 27%  120M 3s
    ##  36350K .......... .......... .......... .......... .......... 27% 9.17M 3s
    ##  36400K .......... .......... .......... .......... .......... 27% 93.3M 3s
    ##  36450K .......... .......... .......... .......... .......... 27% 20.0M 3s
    ##  36500K .......... .......... .......... .......... .......... 27% 51.7M 3s
    ##  36550K .......... .......... .......... .......... .......... 27% 21.3M 3s
    ##  36600K .......... .......... .......... .......... .......... 27% 30.4M 3s
    ##  36650K .......... .......... .......... .......... .......... 27% 17.7M 3s
    ##  36700K .......... .......... .......... .......... .......... 27% 96.4M 3s
    ##  36750K .......... .......... .......... .......... .......... 27% 19.4M 3s
    ##  36800K .......... .......... .......... .......... .......... 27% 91.0M 3s
    ##  36850K .......... .......... .......... .......... .......... 27% 44.1M 3s
    ##  36900K .......... .......... .......... .......... .......... 27% 19.3M 3s
    ##  36950K .......... .......... .......... .......... .......... 27% 38.7M 3s
    ##  37000K .......... .......... .......... .......... .......... 27% 28.4M 3s
    ##  37050K .......... .......... .......... .......... .......... 27% 81.2M 3s
    ##  37100K .......... .......... .......... .......... .......... 27% 28.8M 3s
    ##  37150K .......... .......... .......... .......... .......... 27% 25.8M 3s
    ##  37200K .......... .......... .......... .......... .......... 27% 32.0M 3s
    ##  37250K .......... .......... .......... .......... .......... 27% 31.4M 3s
    ##  37300K .......... .......... .......... .......... .......... 27% 38.0M 3s
    ##  37350K .......... .......... .......... .......... .......... 27% 53.7M 3s
    ##  37400K .......... .......... .......... .......... .......... 27% 13.0M 3s
    ##  37450K .......... .......... .......... .......... .......... 27% 89.6M 3s
    ##  37500K .......... .......... .......... .......... .......... 28%  122M 3s
    ##  37550K .......... .......... .......... .......... .......... 28% 14.2M 3s
    ##  37600K .......... .......... .......... .......... .......... 28%  116M 3s
    ##  37650K .......... .......... .......... .......... .......... 28% 17.4M 3s
    ##  37700K .......... .......... .......... .......... .......... 28% 81.9M 3s
    ##  37750K .......... .......... .......... .......... .......... 28% 49.0M 3s
    ##  37800K .......... .......... .......... .......... .......... 28% 16.5M 3s
    ##  37850K .......... .......... .......... .......... .......... 28% 45.6M 3s
    ##  37900K .......... .......... .......... .......... .......... 28%  106M 3s
    ##  37950K .......... .......... .......... .......... .......... 28% 18.3M 3s
    ##  38000K .......... .......... .......... .......... .......... 28% 42.0M 3s
    ##  38050K .......... .......... .......... .......... .......... 28% 23.0M 3s
    ##  38100K .......... .......... .......... .......... .......... 28% 32.1M 3s
    ##  38150K .......... .......... .......... .......... .......... 28% 85.7M 3s
    ##  38200K .......... .......... .......... .......... .......... 28% 27.5M 3s
    ##  38250K .......... .......... .......... .......... .......... 28% 26.6M 3s
    ##  38300K .......... .......... .......... .......... .......... 28% 46.9M 3s
    ##  38350K .......... .......... .......... .......... .......... 28% 19.7M 3s
    ##  38400K .......... .......... .......... .......... .......... 28% 32.4M 3s
    ##  38450K .......... .......... .......... .......... .......... 28% 83.2M 3s
    ##  38500K .......... .......... .......... .......... .......... 28% 28.5M 3s
    ##  38550K .......... .......... .......... .......... .......... 28% 23.4M 3s
    ##  38600K .......... .......... .......... .......... .......... 28% 26.0M 3s
    ##  38650K .......... .......... .......... .......... .......... 28% 72.0M 3s
    ##  38700K .......... .......... .......... .......... .......... 28% 34.0M 3s
    ##  38750K .......... .......... .......... .......... .......... 28% 22.9M 3s
    ##  38800K .......... .......... .......... .......... .......... 28% 30.8M 3s
    ##  38850K .......... .......... .......... .......... .......... 29% 23.2M 3s
    ##  38900K .......... .......... .......... .......... .......... 29% 94.3M 3s
    ##  38950K .......... .......... .......... .......... .......... 29% 26.8M 3s
    ##  39000K .......... .......... .......... .......... .......... 29% 12.0M 3s
    ##  39050K .......... .......... .......... .......... .......... 29% 79.2M 3s
    ##  39100K .......... .......... .......... .......... .......... 29%  115M 3s
    ##  39150K .......... .......... .......... .......... .......... 29% 13.1M 3s
    ##  39200K .......... .......... .......... .......... .......... 29%  102M 3s
    ##  39250K .......... .......... .......... .......... .......... 29% 15.5M 3s
    ##  39300K .......... .......... .......... .......... .......... 29% 88.0M 3s
    ##  39350K .......... .......... .......... .......... .......... 29% 98.8M 3s
    ##  39400K .......... .......... .......... .......... .......... 29% 14.2M 3s
    ##  39450K .......... .......... .......... .......... .......... 29%  102M 3s
    ##  39500K .......... .......... .......... .......... .......... 29% 14.4M 3s
    ##  39550K .......... .......... .......... .......... .......... 29% 74.0M 3s
    ##  39600K .......... .......... .......... .......... .......... 29% 90.8M 3s
    ##  39650K .......... .......... .......... .......... .......... 29% 18.4M 3s
    ##  39700K .......... .......... .......... .......... .......... 29% 48.7M 3s
    ##  39750K .......... .......... .......... .......... .......... 29% 28.3M 3s
    ##  39800K .......... .......... .......... .......... .......... 29% 48.9M 3s
    ##  39850K .......... .......... .......... .......... .......... 29% 39.9M 3s
    ##  39900K .......... .......... .......... .......... .......... 29% 17.2M 3s
    ##  39950K .......... .......... .......... .......... .......... 29% 38.5M 3s
    ##  40000K .......... .......... .......... .......... .......... 29% 12.5M 3s
    ##  40050K .......... .......... .......... .......... .......... 29% 70.8M 3s
    ##  40100K .......... .......... .......... .......... .......... 29% 98.3M 3s
    ##  40150K .......... .......... .......... .......... .......... 29% 17.8M 3s
    ##  40200K .......... .......... .......... .......... .......... 30% 94.4M 3s
    ##  40250K .......... .......... .......... .......... .......... 30%  105M 3s
    ##  40300K .......... .......... .......... .......... .......... 30% 13.6M 3s
    ##  40350K .......... .......... .......... .......... .......... 30% 96.7M 3s
    ##  40400K .......... .......... .......... .......... .......... 30% 11.4M 3s
    ##  40450K .......... .......... .......... .......... .......... 30% 96.4M 3s
    ##  40500K .......... .......... .......... .......... .......... 30% 12.8M 3s
    ##  40550K .......... .......... .......... .......... .......... 30% 64.7M 3s
    ##  40600K .......... .......... .......... .......... .......... 30% 98.4M 3s
    ##  40650K .......... .......... .......... .......... .......... 30% 17.1M 3s
    ##  40700K .......... .......... .......... .......... .......... 30% 70.9M 3s
    ##  40750K .......... .......... .......... .......... .......... 30% 17.8M 3s
    ##  40800K .......... .......... .......... .......... .......... 30% 79.4M 3s
    ##  40850K .......... .......... .......... .......... .......... 30% 95.3M 3s
    ##  40900K .......... .......... .......... .......... .......... 30% 19.9M 3s
    ##  40950K .......... .......... .......... .......... .......... 30% 62.2M 3s
    ##  41000K .......... .......... .......... .......... .......... 30% 90.7M 3s
    ##  41050K .......... .......... .......... .......... .......... 30% 18.3M 3s
    ##  41100K .......... .......... .......... .......... .......... 30% 29.8M 3s
    ##  41150K .......... .......... .......... .......... .......... 30% 25.8M 3s
    ##  41200K .......... .......... .......... .......... .......... 30% 28.8M 3s
    ##  41250K .......... .......... .......... .......... .......... 30%  101M 3s
    ##  41300K .......... .......... .......... .......... .......... 30% 27.4M 3s
    ##  41350K .......... .......... .......... .......... .......... 30% 26.7M 3s
    ##  41400K .......... .......... .......... .......... .......... 30% 75.2M 3s
    ##  41450K .......... .......... .......... .......... .......... 30% 49.9M 3s
    ##  41500K .......... .......... .......... .......... .......... 30% 18.1M 3s
    ##  41550K .......... .......... .......... .......... .......... 31% 59.1M 3s
    ##  41600K .......... .......... .......... .......... .......... 31% 13.1M 3s
    ##  41650K .......... .......... .......... .......... .......... 31% 99.5M 3s
    ##  41700K .......... .......... .......... .......... .......... 31%  110M 3s
    ##  41750K .......... .......... .......... .......... .......... 31% 14.4M 3s
    ##  41800K .......... .......... .......... .......... .......... 31%  106M 3s
    ##  41850K .......... .......... .......... .......... .......... 31% 15.7M 3s
    ##  41900K .......... .......... .......... .......... .......... 31% 85.1M 3s
    ##  41950K .......... .......... .......... .......... .......... 31% 96.0M 3s
    ##  42000K .......... .......... .......... .......... .......... 31% 20.5M 3s
    ##  42050K .......... .......... .......... .......... .......... 31% 88.4M 3s
    ##  42100K .......... .......... .......... .......... .......... 31% 20.1M 3s
    ##  42150K .......... .......... .......... .......... .......... 31% 89.2M 3s
    ##  42200K .......... .......... .......... .......... .......... 31% 55.2M 3s
    ##  42250K .......... .......... .......... .......... .......... 31% 21.5M 3s
    ##  42300K .......... .......... .......... .......... .......... 31% 58.3M 3s
    ##  42350K .......... .......... .......... .......... .......... 31% 20.5M 3s
    ##  42400K .......... .......... .......... .......... .......... 31% 96.9M 3s
    ##  42450K .......... .......... .......... .......... .......... 31% 45.0M 3s
    ##  42500K .......... .......... .......... .......... .......... 31% 21.2M 3s
    ##  42550K .......... .......... .......... .......... .......... 31% 42.5M 3s
    ##  42600K .......... .......... .......... .......... .......... 31% 67.7M 3s
    ##  42650K .......... .......... .......... .......... .......... 31% 98.3M 3s
    ##  42700K .......... .......... .......... .......... .......... 31% 11.6M 3s
    ##  42750K .......... .......... .......... .......... .......... 31%  101M 3s
    ##  42800K .......... .......... .......... .......... .......... 31%  118M 3s
    ##  42850K .......... .......... .......... .......... .......... 31% 15.3M 3s
    ##  42900K .......... .......... .......... .......... .......... 32%  106M 3s
    ##  42950K .......... .......... .......... .......... .......... 32% 65.0M 3s
    ##  43000K .......... .......... .......... .......... .......... 32% 22.8M 3s
    ##  43050K .......... .......... .......... .......... .......... 32% 28.4M 3s
    ##  43100K .......... .......... .......... .......... .......... 32%  124M 3s
    ##  43150K .......... .......... .......... .......... .......... 32% 38.7M 3s
    ##  43200K .......... .......... .......... .......... .......... 32% 14.9M 3s
    ##  43250K .......... .......... .......... .......... .......... 32%  118M 3s
    ##  43300K .......... .......... .......... .......... .......... 32% 13.9M 3s
    ##  43350K .......... .......... .......... .......... .......... 32% 76.7M 3s
    ##  43400K .......... .......... .......... .......... .......... 32%  103M 3s
    ##  43450K .......... .......... .......... .......... .......... 32% 15.0M 3s
    ##  43500K .......... .......... .......... .......... .......... 32% 87.3M 3s
    ##  43550K .......... .......... .......... .......... .......... 32% 16.8M 3s
    ##  43600K .......... .......... .......... .......... .......... 32% 93.8M 3s
    ##  43650K .......... .......... .......... .......... .......... 32% 56.2M 3s
    ##  43700K .......... .......... .......... .......... .......... 32% 20.5M 3s
    ##  43750K .......... .......... .......... .......... .......... 32% 33.6M 3s
    ##  43800K .......... .......... .......... .......... .......... 32%  144M 3s
    ##  43850K .......... .......... .......... .......... .......... 32% 19.5M 3s
    ##  43900K .......... .......... .......... .......... .......... 32% 39.0M 3s
    ##  43950K .......... .......... .......... .......... .......... 32% 18.8M 3s
    ##  44000K .......... .......... .......... .......... .......... 32% 39.3M 3s
    ##  44050K .......... .......... .......... .......... .......... 32%  135M 3s
    ##  44100K .......... .......... .......... .......... .......... 32% 19.8M 3s
    ##  44150K .......... .......... .......... .......... .......... 32% 38.1M 3s
    ##  44200K .......... .......... .......... .......... .......... 33% 20.1M 3s
    ##  44250K .......... .......... .......... .......... .......... 33%  148M 3s
    ##  44300K .......... .......... .......... .......... .......... 33% 33.4M 3s
    ##  44350K .......... .......... .......... .......... .......... 33% 24.0M 3s
    ##  44400K .......... .......... .......... .......... .......... 33% 24.5M 3s
    ##  44450K .......... .......... .......... .......... .......... 33% 28.7M 3s
    ##  44500K .......... .......... .......... .......... .......... 33%  129M 3s
    ##  44550K .......... .......... .......... .......... .......... 33% 22.2M 3s
    ##  44600K .......... .......... .......... .......... .......... 33% 27.9M 3s
    ##  44650K .......... .......... .......... .......... .......... 33% 21.8M 3s
    ##  44700K .......... .......... .......... .......... .......... 33% 42.5M 3s
    ##  44750K .......... .......... .......... .......... .......... 33% 17.3M 3s
    ##  44800K .......... .......... .......... .......... .......... 33%  117M 3s
    ##  44850K .......... .......... .......... .......... .......... 33% 53.8M 3s
    ##  44900K .......... .......... .......... .......... .......... 33% 13.3M 3s
    ##  44950K .......... .......... .......... .......... .......... 33% 57.8M 3s
    ##  45000K .......... .......... .......... .......... .......... 33%  110M 3s
    ##  45050K .......... .......... .......... .......... .......... 33% 20.2M 3s
    ##  45100K .......... .......... .......... .......... .......... 33% 73.0M 3s
    ##  45150K .......... .......... .......... .......... .......... 33% 16.6M 3s
    ##  45200K .......... .......... .......... .......... .......... 33%  108M 3s
    ##  45250K .......... .......... .......... .......... .......... 33% 99.1M 3s
    ##  45300K .......... .......... .......... .......... .......... 33% 14.3M 3s
    ##  45350K .......... .......... .......... .......... .......... 33% 80.5M 3s
    ##  45400K .......... .......... .......... .......... .......... 33% 14.7M 3s
    ##  45450K .......... .......... .......... .......... .......... 33% 72.1M 3s
    ##  45500K .......... .......... .......... .......... .......... 33% 99.1M 3s
    ##  45550K .......... .......... .......... .......... .......... 34% 15.4M 3s
    ##  45600K .......... .......... .......... .......... .......... 34% 78.7M 3s
    ##  45650K .......... .......... .......... .......... .......... 34% 14.8M 3s
    ##  45700K .......... .......... .......... .......... .......... 34% 69.5M 3s
    ##  45750K .......... .......... .......... .......... .......... 34% 90.2M 3s
    ##  45800K .......... .......... .......... .......... .......... 34% 14.9M 3s
    ##  45850K .......... .......... .......... .......... .......... 34% 95.8M 3s
    ##  45900K .......... .......... .......... .......... .......... 34%  112M 3s
    ##  45950K .......... .......... .......... .......... .......... 34% 14.7M 3s
    ##  46000K .......... .......... .......... .......... .......... 34%  114M 3s
    ##  46050K .......... .......... .......... .......... .......... 34% 14.1M 3s
    ##  46100K .......... .......... .......... .......... .......... 34% 80.9M 3s
    ##  46150K .......... .......... .......... .......... .......... 34% 95.6M 3s
    ##  46200K .......... .......... .......... .......... .......... 34% 19.2M 3s
    ##  46250K .......... .......... .......... .......... .......... 34% 71.4M 3s
    ##  46300K .......... .......... .......... .......... .......... 34% 18.6M 3s
    ##  46350K .......... .......... .......... .......... .......... 34% 26.0M 3s
    ##  46400K .......... .......... .......... .......... .......... 34%  114M 3s
    ##  46450K .......... .......... .......... .......... .......... 34% 17.6M 3s
    ##  46500K .......... .......... .......... .......... .......... 34% 48.0M 3s
    ##  46550K .......... .......... .......... .......... .......... 34% 25.0M 3s
    ##  46600K .......... .......... .......... .......... .......... 34% 87.1M 3s
    ##  46650K .......... .......... .......... .......... .......... 34% 57.6M 3s
    ##  46700K .......... .......... .......... .......... .......... 34% 21.3M 3s
    ##  46750K .......... .......... .......... .......... .......... 34% 81.8M 3s
    ##  46800K .......... .......... .......... .......... .......... 34%  104M 3s
    ##  46850K .......... .......... .......... .......... .......... 34% 17.6M 3s
    ##  46900K .......... .......... .......... .......... .......... 35% 86.6M 3s
    ##  46950K .......... .......... .......... .......... .......... 35%  107M 3s
    ##  47000K .......... .......... .......... .......... .......... 35% 21.5M 3s
    ##  47050K .......... .......... .......... .......... .......... 35% 65.8M 3s
    ##  47100K .......... .......... .......... .......... .......... 35% 80.0M 3s
    ##  47150K .......... .......... .......... .......... .......... 35% 20.3M 3s
    ##  47200K .......... .......... .......... .......... .......... 35% 39.6M 3s
    ##  47250K .......... .......... .......... .......... .......... 35% 27.4M 3s
    ##  47300K .......... .......... .......... .......... .......... 35% 27.3M 3s
    ##  47350K .......... .......... .......... .......... .......... 35% 78.5M 3s
    ##  47400K .......... .......... .......... .......... .......... 35% 37.1M 3s
    ##  47450K .......... .......... .......... .......... .......... 35% 27.3M 3s
    ##  47500K .......... .......... .......... .......... .......... 35% 49.9M 3s
    ##  47550K .......... .......... .......... .......... .......... 35% 20.5M 3s
    ##  47600K .......... .......... .......... .......... .......... 35% 92.3M 3s
    ##  47650K .......... .......... .......... .......... .......... 35% 44.2M 3s
    ##  47700K .......... .......... .......... .......... .......... 35% 22.2M 3s
    ##  47750K .......... .......... .......... .......... .......... 35% 44.0M 3s
    ##  47800K .......... .......... .......... .......... .......... 35% 88.7M 3s
    ##  47850K .......... .......... .......... .......... .......... 35% 27.1M 3s
    ##  47900K .......... .......... .......... .......... .......... 35% 32.2M 3s
    ##  47950K .......... .......... .......... .......... .......... 35% 26.6M 3s
    ##  48000K .......... .......... .......... .......... .......... 35% 26.0M 3s
    ##  48050K .......... .......... .......... .......... .......... 35% 98.6M 3s
    ##  48100K .......... .......... .......... .......... .......... 35% 55.2M 3s
    ##  48150K .......... .......... .......... .......... .......... 35% 17.6M 3s
    ##  48200K .......... .......... .......... .......... .......... 35% 65.7M 3s
    ##  48250K .......... .......... .......... .......... .......... 36%  114M 3s
    ##  48300K .......... .......... .......... .......... .......... 36% 14.8M 3s
    ##  48350K .......... .......... .......... .......... .......... 36% 73.3M 3s
    ##  48400K .......... .......... .......... .......... .......... 36% 18.7M 3s
    ##  48450K .......... .......... .......... .......... .......... 36% 83.3M 3s
    ##  48500K .......... .......... .......... .......... .......... 36%  103M 3s
    ##  48550K .......... .......... .......... .......... .......... 36% 16.4M 3s
    ##  48600K .......... .......... .......... .......... .......... 36%  112M 3s
    ##  48650K .......... .......... .......... .......... .......... 36% 19.3M 3s
    ##  48700K .......... .......... .......... .......... .......... 36% 85.7M 3s
    ##  48750K .......... .......... .......... .......... .......... 36% 67.7M 3s
    ##  48800K .......... .......... .......... .......... .......... 36% 20.2M 3s
    ##  48850K .......... .......... .......... .......... .......... 36% 42.3M 3s
    ##  48900K .......... .......... .......... .......... .......... 36% 22.5M 3s
    ##  48950K .......... .......... .......... .......... .......... 36% 41.0M 3s
    ##  49000K .......... .......... .......... .......... .......... 36%  109M 3s
    ##  49050K .......... .......... .......... .......... .......... 36% 20.1M 3s
    ##  49100K .......... .......... .......... .......... .......... 36% 99.6M 3s
    ##  49150K .......... .......... .......... .......... .......... 36% 18.3M 3s
    ##  49200K .......... .......... .......... .......... .......... 36% 97.8M 3s
    ##  49250K .......... .......... .......... .......... .......... 36% 81.3M 3s
    ##  49300K .......... .......... .......... .......... .......... 36% 21.4M 3s
    ##  49350K .......... .......... .......... .......... .......... 36% 48.5M 3s
    ##  49400K .......... .......... .......... .......... .......... 36% 28.4M 3s
    ##  49450K .......... .......... .......... .......... .......... 36% 85.4M 3s
    ##  49500K .......... .......... .......... .......... .......... 36% 27.8M 3s
    ##  49550K .......... .......... .......... .......... .......... 36% 38.0M 3s
    ##  49600K .......... .......... .......... .......... .......... 37% 22.9M 3s
    ##  49650K .......... .......... .......... .......... .......... 37% 58.9M 3s
    ##  49700K .......... .......... .......... .......... .......... 37% 95.6M 3s
    ##  49750K .......... .......... .......... .......... .......... 37% 23.8M 3s
    ##  49800K .......... .......... .......... .......... .......... 37% 7.71M 3s
    ##  49850K .......... .......... .......... .......... .......... 37%  108M 3s
    ##  49900K .......... .......... .......... .......... .......... 37%  124M 3s
    ##  49950K .......... .......... .......... .......... .......... 37% 14.0M 3s
    ##  50000K .......... .......... .......... .......... .......... 37% 85.8M 3s
    ##  50050K .......... .......... .......... .......... .......... 37%  109M 3s
    ##  50100K .......... .......... .......... .......... .......... 37% 15.1M 3s
    ##  50150K .......... .......... .......... .......... .......... 37% 83.9M 3s
    ##  50200K .......... .......... .......... .......... .......... 37%  120M 3s
    ##  50250K .......... .......... .......... .......... .......... 37% 12.9M 3s
    ##  50300K .......... .......... .......... .......... .......... 37%  116M 3s
    ##  50350K .......... .......... .......... .......... .......... 37% 14.1M 3s
    ##  50400K .......... .......... .......... .......... .......... 37% 82.7M 3s
    ##  50450K .......... .......... .......... .......... .......... 37% 73.3M 3s
    ##  50500K .......... .......... .......... .......... .......... 37%  116M 3s
    ##  50550K .......... .......... .......... .......... .......... 37% 29.2M 3s
    ##  50600K .......... .......... .......... .......... .......... 37% 90.5M 3s
    ##  50650K .......... .......... .......... .......... .......... 37% 67.7M 3s
    ##  50700K .......... .......... .......... .......... .......... 37% 17.0M 3s
    ##  50750K .......... .......... .......... .......... .......... 37% 72.8M 3s
    ##  50800K .......... .......... .......... .......... .......... 37% 99.6M 3s
    ##  50850K .......... .......... .......... .......... .......... 37% 93.3M 3s
    ##  50900K .......... .......... .......... .......... .......... 38% 17.7M 3s
    ##  50950K .......... .......... .......... .......... .......... 38% 81.0M 3s
    ##  51000K .......... .......... .......... .......... .......... 38% 46.5M 3s
    ##  51050K .......... .......... .......... .......... .......... 38% 44.5M 3s
    ##  51100K .......... .......... .......... .......... .......... 38% 73.3M 3s
    ##  51150K .......... .......... .......... .......... .......... 38% 22.4M 3s
    ##  51200K .......... .......... .......... .......... .......... 38% 97.4M 3s
    ##  51250K .......... .......... .......... .......... .......... 38% 40.7M 3s
    ##  51300K .......... .......... .......... .......... .......... 38% 22.3M 3s
    ##  51350K .......... .......... .......... .......... .......... 38% 69.2M 3s
    ##  51400K .......... .......... .......... .......... .......... 38%  114M 3s
    ##  51450K .......... .......... .......... .......... .......... 38% 57.5M 3s
    ##  51500K .......... .......... .......... .......... .......... 38% 30.8M 3s
    ##  51550K .......... .......... .......... .......... .......... 38% 68.0M 3s
    ##  51600K .......... .......... .......... .......... .......... 38% 96.5M 3s
    ##  51650K .......... .......... .......... .......... .......... 38% 48.4M 3s
    ##  51700K .......... .......... .......... .......... .......... 38% 15.9M 3s
    ##  51750K .......... .......... .......... .......... .......... 38% 73.4M 3s
    ##  51800K .......... .......... .......... .......... .......... 38% 97.2M 3s
    ##  51850K .......... .......... .......... .......... .......... 38%  114M 3s
    ##  51900K .......... .......... .......... .......... .......... 38% 14.1M 3s
    ##  51950K .......... .......... .......... .......... .......... 38% 83.4M 3s
    ##  52000K .......... .......... .......... .......... .......... 38% 28.6M 3s
    ##  52050K .......... .......... .......... .......... .......... 38% 70.6M 3s
    ##  52100K .......... .......... .......... .......... .......... 38% 26.2M 3s
    ##  52150K .......... .......... .......... .......... .......... 38%  118M 3s
    ##  52200K .......... .......... .......... .......... .......... 38% 90.9M 3s
    ##  52250K .......... .......... .......... .......... .......... 39%  137M 3s
    ##  52300K .......... .......... .......... .......... .......... 39% 18.9M 3s
    ##  52350K .......... .......... .......... .......... .......... 39% 71.6M 3s
    ##  52400K .......... .......... .......... .......... .......... 39% 68.8M 3s
    ##  52450K .......... .......... .......... .......... .......... 39% 21.2M 3s
    ##  52500K .......... .......... .......... .......... .......... 39% 98.6M 3s
    ##  52550K .......... .......... .......... .......... .......... 39% 84.5M 3s
    ##  52600K .......... .......... .......... .......... .......... 39% 57.6M 3s
    ##  52650K .......... .......... .......... .......... .......... 39% 21.7M 3s
    ##  52700K .......... .......... .......... .......... .......... 39%  133M 3s
    ##  52750K .......... .......... .......... .......... .......... 39% 61.7M 3s
    ##  52800K .......... .......... .......... .......... .......... 39% 58.3M 3s
    ##  52850K .......... .......... .......... .......... .......... 39% 13.7M 3s
    ##  52900K .......... .......... .......... .......... .......... 39%  121M 3s
    ##  52950K .......... .......... .......... .......... .......... 39%  106M 3s
    ##  53000K .......... .......... .......... .......... .......... 39%  129M 3s
    ##  53050K .......... .......... .......... .......... .......... 39% 11.9M 3s
    ##  53100K .......... .......... .......... .......... .......... 39% 71.8M 3s
    ##  53150K .......... .......... .......... .......... .......... 39% 84.3M 3s
    ##  53200K .......... .......... .......... .......... .......... 39% 18.1M 3s
    ##  53250K .......... .......... .......... .......... .......... 39%  102M 3s
    ##  53300K .......... .......... .......... .......... .......... 39% 5.93M 3s
    ##  53350K .......... .......... .......... .......... .......... 39%  119M 3s
    ##  53400K .......... .......... .......... .......... .......... 39%  138M 3s
    ##  53450K .......... .......... .......... .......... .......... 39%  147M 3s
    ##  53500K .......... .......... .......... .......... .......... 39%  126M 3s
    ##  53550K .......... .......... .......... .......... .......... 39%  121M 3s
    ##  53600K .......... .......... .......... .......... .......... 40% 28.3M 3s
    ##  53650K .......... .......... .......... .......... .......... 40% 5.45M 3s
    ##  53700K .......... .......... .......... .......... .......... 40%  131M 3s
    ##  53750K .......... .......... .......... .......... .......... 40%  113M 3s
    ##  53800K .......... .......... .......... .......... .......... 40%  156M 3s
    ##  53850K .......... .......... .......... .......... .......... 40%  148M 3s
    ##  53900K .......... .......... .......... .......... .......... 40%  147M 3s
    ##  53950K .......... .......... .......... .......... .......... 40% 26.4M 3s
    ##  54000K .......... .......... .......... .......... .......... 40% 80.1M 2s
    ##  54050K .......... .......... .......... .......... .......... 40% 16.7M 2s
    ##  54100K .......... .......... .......... .......... .......... 40% 17.8M 2s
    ##  54150K .......... .......... .......... .......... .......... 40% 64.3M 2s
    ##  54200K .......... .......... .......... .......... .......... 40% 19.2M 2s
    ##  54250K .......... .......... .......... .......... .......... 40% 69.0M 2s
    ##  54300K .......... .......... .......... .......... .......... 40% 18.7M 2s
    ##  54350K .......... .......... .......... .......... .......... 40% 24.6M 2s
    ##  54400K .......... .......... .......... .......... .......... 40% 13.5M 2s
    ##  54450K .......... .......... .......... .......... .......... 40%  156M 2s
    ##  54500K .......... .......... .......... .......... .......... 40% 11.2M 2s
    ##  54550K .......... .......... .......... .......... .......... 40%  148M 2s
    ##  54600K .......... .......... .......... .......... .......... 40% 13.2M 2s
    ##  54650K .......... .......... .......... .......... .......... 40%  145M 2s
    ##  54700K .......... .......... .......... .......... .......... 40% 13.2M 2s
    ##  54750K .......... .......... .......... .......... .......... 40%  107M 2s
    ##  54800K .......... .......... .......... .......... .......... 40% 11.9M 2s
    ##  54850K .......... .......... .......... .......... .......... 40%  164M 2s
    ##  54900K .......... .......... .......... .......... .......... 40% 13.5M 2s
    ##  54950K .......... .......... .......... .......... .......... 41%  102M 2s
    ##  55000K .......... .......... .......... .......... .......... 41% 17.6M 2s
    ##  55050K .......... .......... .......... .......... .......... 41% 63.6M 2s
    ##  55100K .......... .......... .......... .......... .......... 41% 12.8M 2s
    ##  55150K .......... .......... .......... .......... .......... 41% 48.1M 2s
    ##  55200K .......... .......... .......... .......... .......... 41% 12.3M 2s
    ##  55250K .......... .......... .......... .......... .......... 41% 87.6M 2s
    ##  55300K .......... .......... .......... .......... .......... 41% 12.4M 2s
    ##  55350K .......... .......... .......... .......... .......... 41% 93.9M 2s
    ##  55400K .......... .......... .......... .......... .......... 41% 12.1M 2s
    ##  55450K .......... .......... .......... .......... .......... 41%  102M 2s
    ##  55500K .......... .......... .......... .......... .......... 41% 14.4M 2s
    ##  55550K .......... .......... .......... .......... .......... 41% 81.7M 2s
    ##  55600K .......... .......... .......... .......... .......... 41% 13.5M 2s
    ##  55650K .......... .......... .......... .......... .......... 41%  103M 2s
    ##  55700K .......... .......... .......... .......... .......... 41% 12.4M 2s
    ##  55750K .......... .......... .......... .......... .......... 41% 96.3M 2s
    ##  55800K .......... .......... .......... .......... .......... 41% 13.4M 2s
    ##  55850K .......... .......... .......... .......... .......... 41%  105M 2s
    ##  55900K .......... .......... .......... .......... .......... 41% 12.2M 2s
    ##  55950K .......... .......... .......... .......... .......... 41% 87.7M 2s
    ##  56000K .......... .......... .......... .......... .......... 41% 16.9M 2s
    ##  56050K .......... .......... .......... .......... .......... 41% 94.5M 2s
    ##  56100K .......... .......... .......... .......... .......... 41% 14.7M 2s
    ##  56150K .......... .......... .......... .......... .......... 41% 92.0M 2s
    ##  56200K .......... .......... .......... .......... .......... 41% 15.6M 2s
    ##  56250K .......... .......... .......... .......... .......... 41% 82.9M 2s
    ##  56300K .......... .......... .......... .......... .......... 42% 93.8M 2s
    ##  56350K .......... .......... .......... .......... .......... 42% 13.4M 2s
    ##  56400K .......... .......... .......... .......... .......... 42% 15.3M 2s
    ##  56450K .......... .......... .......... .......... .......... 42% 95.1M 2s
    ##  56500K .......... .......... .......... .......... .......... 42% 18.5M 2s
    ##  56550K .......... .......... .......... .......... .......... 42% 46.4M 2s
    ##  56600K .......... .......... .......... .......... .......... 42% 12.4M 2s
    ##  56650K .......... .......... .......... .......... .......... 42% 86.3M 2s
    ##  56700K .......... .......... .......... .......... .......... 42% 18.4M 2s
    ##  56750K .......... .......... .......... .......... .......... 42% 59.4M 2s
    ##  56800K .......... .......... .......... .......... .......... 42% 85.1M 2s
    ##  56850K .......... .......... .......... .......... .......... 42% 15.0M 2s
    ##  56900K .......... .......... .......... .......... .......... 42% 23.0M 2s
    ##  56950K .......... .......... .......... .......... .......... 42% 33.9M 2s
    ##  57000K .......... .......... .......... .......... .......... 42% 18.3M 2s
    ##  57050K .......... .......... .......... .......... .......... 42%  104M 2s
    ##  57100K .......... .......... .......... .......... .......... 42% 32.1M 2s
    ##  57150K .......... .......... .......... .......... .......... 42% 20.6M 2s
    ##  57200K .......... .......... .......... .......... .......... 42% 59.7M 2s
    ##  57250K .......... .......... .......... .......... .......... 42% 21.5M 2s
    ##  57300K .......... .......... .......... .......... .......... 42% 42.8M 2s
    ##  57350K .......... .......... .......... .......... .......... 42% 16.5M 2s
    ##  57400K .......... .......... .......... .......... .......... 42% 12.3M 2s
    ##  57450K .......... .......... .......... .......... .......... 42% 88.4M 2s
    ##  57500K .......... .......... .......... .......... .......... 42%  109M 2s
    ##  57550K .......... .......... .......... .......... .......... 42% 14.5M 2s
    ##  57600K .......... .......... .......... .......... .......... 43% 11.4M 2s
    ##  57650K .......... .......... .......... .......... .......... 43%  115M 2s
    ##  57700K .......... .......... .......... .......... .......... 43%  121M 2s
    ##  57750K .......... .......... .......... .......... .......... 43% 12.2M 2s
    ##  57800K .......... .......... .......... .......... .......... 43% 97.4M 2s
    ##  57850K .......... .......... .......... .......... .......... 43% 13.9M 2s
    ##  57900K .......... .......... .......... .......... .......... 43% 84.9M 2s
    ##  57950K .......... .......... .......... .......... .......... 43% 20.1M 2s
    ##  58000K .......... .......... .......... .......... .......... 43% 85.9M 2s
    ##  58050K .......... .......... .......... .......... .......... 43% 94.7M 2s
    ##  58100K .......... .......... .......... .......... .......... 43% 13.0M 2s
    ##  58150K .......... .......... .......... .......... .......... 43% 80.9M 2s
    ##  58200K .......... .......... .......... .......... .......... 43%  118M 2s
    ##  58250K .......... .......... .......... .......... .......... 43% 16.3M 2s
    ##  58300K .......... .......... .......... .......... .......... 43% 13.0M 2s
    ##  58350K .......... .......... .......... .......... .......... 43% 59.0M 2s
    ##  58400K .......... .......... .......... .......... .......... 43%  115M 2s
    ##  58450K .......... .......... .......... .......... .......... 43% 12.8M 2s
    ##  58500K .......... .......... .......... .......... .......... 43%  108M 2s
    ##  58550K .......... .......... .......... .......... .......... 43% 3.10M 2s
    ##  58600K .......... .......... .......... .......... .......... 43%  108M 2s
    ##  58650K .......... .......... .......... .......... .......... 43% 69.0M 2s
    ##  58700K .......... .......... .......... .......... .......... 43% 14.6M 2s
    ##  58750K .......... .......... .......... .......... .......... 43% 91.9M 2s
    ##  58800K .......... .......... .......... .......... .......... 43%  105M 2s
    ##  58850K .......... .......... .......... .......... .......... 43%  106M 2s
    ##  58900K .......... .......... .......... .......... .......... 43% 17.9M 2s
    ##  58950K .......... .......... .......... .......... .......... 44% 92.3M 2s
    ##  59000K .......... .......... .......... .......... .......... 44% 13.1M 2s
    ##  59050K .......... .......... .......... .......... .......... 44% 84.1M 2s
    ##  59100K .......... .......... .......... .......... .......... 44% 12.7M 2s
    ##  59150K .......... .......... .......... .......... .......... 44% 11.8M 2s
    ##  59200K .......... .......... .......... .......... .......... 44% 13.1M 2s
    ##  59250K .......... .......... .......... .......... .......... 44% 13.2M 2s
    ##  59300K .......... .......... .......... .......... .......... 44% 71.9M 2s
    ##  59350K .......... .......... .......... .......... .......... 44% 15.5M 2s
    ##  59400K .......... .......... .......... .......... .......... 44% 16.0M 2s
    ##  59450K .......... .......... .......... .......... .......... 44% 13.1M 2s
    ##  59500K .......... .......... .......... .......... .......... 44% 13.1M 2s
    ##  59550K .......... .......... .......... .......... .......... 44% 12.9M 2s
    ##  59600K .......... .......... .......... .......... .......... 44% 85.0M 2s
    ##  59650K .......... .......... .......... .......... .......... 44% 12.3M 2s
    ##  59700K .......... .......... .......... .......... .......... 44% 12.5M 2s
    ##  59750K .......... .......... .......... .......... .......... 44% 11.5M 2s
    ##  59800K .......... .......... .......... .......... .......... 44% 88.6M 2s
    ##  59850K .......... .......... .......... .......... .......... 44% 12.9M 2s
    ##  59900K .......... .......... .......... .......... .......... 44% 15.3M 2s
    ##  59950K .......... .......... .......... .......... .......... 44% 69.1M 2s
    ##  60000K .......... .......... .......... .......... .......... 44% 16.2M 2s
    ##  60050K .......... .......... .......... .......... .......... 44% 77.0M 2s
    ##  60100K .......... .......... .......... .......... .......... 44% 16.3M 2s
    ##  60150K .......... .......... .......... .......... .......... 44% 62.9M 2s
    ##  60200K .......... .......... .......... .......... .......... 44% 11.8M 2s
    ##  60250K .......... .......... .......... .......... .......... 44%  102M 2s
    ##  60300K .......... .......... .......... .......... .......... 45% 12.7M 2s
    ##  60350K .......... .......... .......... .......... .......... 45% 12.4M 2s
    ##  60400K .......... .......... .......... .......... .......... 45% 85.3M 2s
    ##  60450K .......... .......... .......... .......... .......... 45% 12.0M 2s
    ##  60500K .......... .......... .......... .......... .......... 45% 93.1M 2s
    ##  60550K .......... .......... .......... .......... .......... 45% 12.4M 2s
    ##  60600K .......... .......... .......... .......... .......... 45% 96.6M 2s
    ##  60650K .......... .......... .......... .......... .......... 45% 17.5M 2s
    ##  60700K .......... .......... .......... .......... .......... 45% 56.4M 2s
    ##  60750K .......... .......... .......... .......... .......... 45% 11.1M 2s
    ##  60800K .......... .......... .......... .......... .......... 45%  106M 2s
    ##  60850K .......... .......... .......... .......... .......... 45% 11.5M 2s
    ##  60900K .......... .......... .......... .......... .......... 45%  101M 2s
    ##  60950K .......... .......... .......... .......... .......... 45% 12.0M 2s
    ##  61000K .......... .......... .......... .......... .......... 45% 95.9M 2s
    ##  61050K .......... .......... .......... .......... .......... 45% 11.9M 2s
    ##  61100K .......... .......... .......... .......... .......... 45%  109M 2s
    ##  61150K .......... .......... .......... .......... .......... 45% 12.5M 2s
    ##  61200K .......... .......... .......... .......... .......... 45%  102M 2s
    ##  61250K .......... .......... .......... .......... .......... 45% 11.4M 2s
    ##  61300K .......... .......... .......... .......... .......... 45%  108M 2s
    ##  61350K .......... .......... .......... .......... .......... 45% 12.4M 2s
    ##  61400K .......... .......... .......... .......... .......... 45%  101M 2s
    ##  61450K .......... .......... .......... .......... .......... 45% 11.6M 2s
    ##  61500K .......... .......... .......... .......... .......... 45%  105M 2s
    ##  61550K .......... .......... .......... .......... .......... 45% 12.6M 2s
    ##  61600K .......... .......... .......... .......... .......... 45% 14.3M 2s
    ##  61650K .......... .......... .......... .......... .......... 46% 72.3M 2s
    ##  61700K .......... .......... .......... .......... .......... 46% 19.2M 2s
    ##  61750K .......... .......... .......... .......... .......... 46% 53.2M 2s
    ##  61800K .......... .......... .......... .......... .......... 46% 15.6M 2s
    ##  61850K .......... .......... .......... .......... .......... 46% 91.7M 2s
    ##  61900K .......... .......... .......... .......... .......... 46% 23.7M 2s
    ##  61950K .......... .......... .......... .......... .......... 46% 26.7M 2s
    ##  62000K .......... .......... .......... .......... .......... 46% 12.0M 2s
    ##  62050K .......... .......... .......... .......... .......... 46% 98.9M 2s
    ##  62100K .......... .......... .......... .......... .......... 46% 15.3M 2s
    ##  62150K .......... .......... .......... .......... .......... 46% 68.6M 2s
    ##  62200K .......... .......... .......... .......... .......... 46% 75.6M 2s
    ##  62250K .......... .......... .......... .......... .......... 46% 19.7M 2s
    ##  62300K .......... .......... .......... .......... .......... 46% 13.9M 2s
    ##  62350K .......... .......... .......... .......... .......... 46% 73.1M 2s
    ##  62400K .......... .......... .......... .......... .......... 46% 13.8M 2s
    ##  62450K .......... .......... .......... .......... .......... 46%  106M 2s
    ##  62500K .......... .......... .......... .......... .......... 46%  115M 2s
    ##  62550K .......... .......... .......... .......... .......... 46% 17.8M 2s
    ##  62600K .......... .......... .......... .......... .......... 46% 59.1M 2s
    ##  62650K .......... .......... .......... .......... .......... 46% 20.3M 2s
    ##  62700K .......... .......... .......... .......... .......... 46% 60.3M 2s
    ##  62750K .......... .......... .......... .......... .......... 46% 13.6M 2s
    ##  62800K .......... .......... .......... .......... .......... 46% 22.0M 2s
    ##  62850K .......... .......... .......... .......... .......... 46% 31.1M 2s
    ##  62900K .......... .......... .......... .......... .......... 46%  105M 2s
    ##  62950K .......... .......... .......... .......... .......... 46% 13.3M 2s
    ##  63000K .......... .......... .......... .......... .......... 47%  102M 2s
    ##  63050K .......... .......... .......... .......... .......... 47% 22.3M 2s
    ##  63100K .......... .......... .......... .......... .......... 47% 30.8M 2s
    ##  63150K .......... .......... .......... .......... .......... 47% 27.4M 2s
    ##  63200K .......... .......... .......... .......... .......... 47% 20.4M 2s
    ##  63250K .......... .......... .......... .......... .......... 47% 11.2M 2s
    ##  63300K .......... .......... .......... .......... .......... 47%  101M 2s
    ##  63350K .......... .......... .......... .......... .......... 47% 12.1M 2s
    ##  63400K .......... .......... .......... .......... .......... 47%  111M 2s
    ##  63450K .......... .......... .......... .......... .......... 47% 11.6M 2s
    ##  63500K .......... .......... .......... .......... .......... 47% 95.8M 2s
    ##  63550K .......... .......... .......... .......... .......... 47% 17.8M 2s
    ##  63600K .......... .......... .......... .......... .......... 47% 90.8M 2s
    ##  63650K .......... .......... .......... .......... .......... 47% 66.5M 2s
    ##  63700K .......... .......... .......... .......... .......... 47% 21.0M 2s
    ##  63750K .......... .......... .......... .......... .......... 47% 49.0M 2s
    ##  63800K .......... .......... .......... .......... .......... 47% 25.4M 2s
    ##  63850K .......... .......... .......... .......... .......... 47% 35.6M 2s
    ##  63900K .......... .......... .......... .......... .......... 47% 20.4M 2s
    ##  63950K .......... .......... .......... .......... .......... 47% 25.1M 2s
    ##  64000K .......... .......... .......... .......... .......... 47% 30.2M 2s
    ##  64050K .......... .......... .......... .......... .......... 47% 20.1M 2s
    ##  64100K .......... .......... .......... .......... .......... 47% 95.3M 2s
    ##  64150K .......... .......... .......... .......... .......... 47% 61.8M 2s
    ##  64200K .......... .......... .......... .......... .......... 47% 14.5M 2s
    ##  64250K .......... .......... .......... .......... .......... 47%  114M 2s
    ##  64300K .......... .......... .......... .......... .......... 47% 10.9M 2s
    ##  64350K .......... .......... .......... .......... .......... 48% 96.7M 2s
    ##  64400K .......... .......... .......... .......... .......... 48% 18.8M 2s
    ##  64450K .......... .......... .......... .......... .......... 48% 63.0M 2s
    ##  64500K .......... .......... .......... .......... .......... 48% 27.3M 2s
    ##  64550K .......... .......... .......... .......... .......... 48% 20.2M 2s
    ##  64600K .......... .......... .......... .......... .......... 48% 40.5M 2s
    ##  64650K .......... .......... .......... .......... .......... 48% 16.8M 2s
    ##  64700K .......... .......... .......... .......... .......... 48%  106M 2s
    ##  64750K .......... .......... .......... .......... .......... 48% 9.99M 2s
    ##  64800K .......... .......... .......... .......... .......... 48% 85.3M 2s
    ##  64850K .......... .......... .......... .......... .......... 48% 84.4M 2s
    ##  64900K .......... .......... .......... .......... .......... 48% 14.6M 2s
    ##  64950K .......... .......... .......... .......... .......... 48% 85.5M 2s
    ##  65000K .......... .......... .......... .......... .......... 48%  106M 2s
    ##  65050K .......... .......... .......... .......... .......... 48% 13.6M 2s
    ##  65100K .......... .......... .......... .......... .......... 48%  108M 2s
    ##  65150K .......... .......... .......... .......... .......... 48% 14.0M 2s
    ##  65200K .......... .......... .......... .......... .......... 48% 84.9M 2s
    ##  65250K .......... .......... .......... .......... .......... 48%  104M 2s
    ##  65300K .......... .......... .......... .......... .......... 48% 14.7M 2s
    ##  65350K .......... .......... .......... .......... .......... 48% 86.2M 2s
    ##  65400K .......... .......... .......... .......... .......... 48%  107M 2s
    ##  65450K .......... .......... .......... .......... .......... 48% 17.2M 2s
    ##  65500K .......... .......... .......... .......... .......... 48% 81.5M 2s
    ##  65550K .......... .......... .......... .......... .......... 48% 17.3M 2s
    ##  65600K .......... .......... .......... .......... .......... 48% 65.6M 2s
    ##  65650K .......... .......... .......... .......... .......... 49% 86.4M 2s
    ##  65700K .......... .......... .......... .......... .......... 49% 19.9M 2s
    ##  65750K .......... .......... .......... .......... .......... 49% 78.0M 2s
    ##  65800K .......... .......... .......... .......... .......... 49% 17.9M 2s
    ##  65850K .......... .......... .......... .......... .......... 49% 31.0M 2s
    ##  65900K .......... .......... .......... .......... .......... 49%  120M 2s
    ##  65950K .......... .......... .......... .......... .......... 49% 23.0M 2s
    ##  66000K .......... .......... .......... .......... .......... 49% 25.0M 2s
    ##  66050K .......... .......... .......... .......... .......... 49% 29.1M 2s
    ##  66100K .......... .......... .......... .......... .......... 49% 22.6M 2s
    ##  66150K .......... .......... .......... .......... .......... 49% 90.3M 2s
    ##  66200K .......... .......... .......... .......... .......... 49% 63.2M 2s
    ##  66250K .......... .......... .......... .......... .......... 49% 18.9M 2s
    ##  66300K .......... .......... .......... .......... .......... 49% 52.7M 2s
    ##  66350K .......... .......... .......... .......... .......... 49% 32.7M 2s
    ##  66400K .......... .......... .......... .......... .......... 49% 29.3M 2s
    ##  66450K .......... .......... .......... .......... .......... 49% 71.7M 2s
    ##  66500K .......... .......... .......... .......... .......... 49% 6.00M 2s
    ##  66550K .......... .......... .......... .......... .......... 49%  104M 2s
    ##  66600K .......... .......... .......... .......... .......... 49%  107M 2s
    ##  66650K .......... .......... .......... .......... .......... 49% 19.6M 2s
    ##  66700K .......... .......... .......... .......... .......... 49% 84.1M 2s
    ##  66750K .......... .......... .......... .......... .......... 49% 16.5M 2s
    ##  66800K .......... .......... .......... .......... .......... 49% 79.3M 2s
    ##  66850K .......... .......... .......... .......... .......... 49% 66.1M 2s
    ##  66900K .......... .......... .......... .......... .......... 49% 15.6M 2s
    ##  66950K .......... .......... .......... .......... .......... 49% 78.8M 2s
    ##  67000K .......... .......... .......... .......... .......... 50% 13.8M 2s
    ##  67050K .......... .......... .......... .......... .......... 50% 69.5M 2s
    ##  67100K .......... .......... .......... .......... .......... 50%  110M 2s
    ##  67150K .......... .......... .......... .......... .......... 50% 14.7M 2s
    ##  67200K .......... .......... .......... .......... .......... 50% 94.2M 2s
    ##  67250K .......... .......... .......... .......... .......... 50% 17.9M 2s
    ##  67300K .......... .......... .......... .......... .......... 50% 91.5M 2s
    ##  67350K .......... .......... .......... .......... .......... 50%  117M 2s
    ##  67400K .......... .......... .......... .......... .......... 50% 16.2M 2s
    ##  67450K .......... .......... .......... .......... .......... 50%  100M 2s
    ##  67500K .......... .......... .......... .......... .......... 50%  133M 2s
    ##  67550K .......... .......... .......... .......... .......... 50% 17.0M 2s
    ##  67600K .......... .......... .......... .......... .......... 50%  109M 2s
    ##  67650K .......... .......... .......... .......... .......... 50% 13.9M 2s
    ##  67700K .......... .......... .......... .......... .......... 50%  121M 2s
    ##  67750K .......... .......... .......... .......... .......... 50%  109M 2s
    ##  67800K .......... .......... .......... .......... .......... 50% 14.8M 2s
    ##  67850K .......... .......... .......... .......... .......... 50%  127M 2s
    ##  67900K .......... .......... .......... .......... .......... 50% 15.1M 2s
    ##  67950K .......... .......... .......... .......... .......... 50% 89.6M 2s
    ##  68000K .......... .......... .......... .......... .......... 50% 18.7M 2s
    ##  68050K .......... .......... .......... .......... .......... 50% 84.9M 2s
    ##  68100K .......... .......... .......... .......... .......... 50%  119M 2s
    ##  68150K .......... .......... .......... .......... .......... 50% 16.1M 2s
    ##  68200K .......... .......... .......... .......... .......... 50%  122M 2s
    ##  68250K .......... .......... .......... .......... .......... 50%  127M 2s
    ##  68300K .......... .......... .......... .......... .......... 50% 18.6M 2s
    ##  68350K .......... .......... .......... .......... .......... 51% 72.5M 2s
    ##  68400K .......... .......... .......... .......... .......... 51% 23.3M 2s
    ##  68450K .......... .......... .......... .......... .......... 51% 34.1M 2s
    ##  68500K .......... .......... .......... .......... .......... 51% 72.6M 2s
    ##  68550K .......... .......... .......... .......... .......... 51% 22.7M 2s
    ##  68600K .......... .......... .......... .......... .......... 51% 49.2M 2s
    ##  68650K .......... .......... .......... .......... .......... 51% 16.9M 2s
    ##  68700K .......... .......... .......... .......... .......... 51%  109M 2s
    ##  68750K .......... .......... .......... .......... .......... 51% 58.1M 2s
    ##  68800K .......... .......... .......... .......... .......... 51% 12.7M 2s
    ##  68850K .......... .......... .......... .......... .......... 51% 81.8M 2s
    ##  68900K .......... .......... .......... .......... .......... 51% 19.4M 2s
    ##  68950K .......... .......... .......... .......... .......... 51% 75.2M 2s
    ##  69000K .......... .......... .......... .......... .......... 51% 51.6M 2s
    ##  69050K .......... .......... .......... .......... .......... 51% 24.8M 2s
    ##  69100K .......... .......... .......... .......... .......... 51% 29.1M 2s
    ##  69150K .......... .......... .......... .......... .......... 51% 26.1M 2s
    ##  69200K .......... .......... .......... .......... .......... 51% 95.3M 2s
    ##  69250K .......... .......... .......... .......... .......... 51% 30.0M 2s
    ##  69300K .......... .......... .......... .......... .......... 51% 23.8M 2s
    ##  69350K .......... .......... .......... .......... .......... 51% 33.2M 2s
    ##  69400K .......... .......... .......... .......... .......... 51%  118M 2s
    ##  69450K .......... .......... .......... .......... .......... 51% 18.8M 2s
    ##  69500K .......... .......... .......... .......... .......... 51% 41.2M 2s
    ##  69550K .......... .......... .......... .......... .......... 51% 13.2M 2s
    ##  69600K .......... .......... .......... .......... .......... 51%  112M 2s
    ##  69650K .......... .......... .......... .......... .......... 51%  123M 2s
    ##  69700K .......... .......... .......... .......... .......... 52% 13.0M 2s
    ##  69750K .......... .......... .......... .......... .......... 52%  106M 2s
    ##  69800K .......... .......... .......... .......... .......... 52% 12.8M 2s
    ##  69850K .......... .......... .......... .......... .......... 52%  101M 2s
    ##  69900K .......... .......... .......... .......... .......... 52%  121M 2s
    ##  69950K .......... .......... .......... .......... .......... 52% 13.2M 2s
    ##  70000K .......... .......... .......... .......... .......... 52%  126M 2s
    ##  70050K .......... .......... .......... .......... .......... 52% 14.2M 2s
    ##  70100K .......... .......... .......... .......... .......... 52% 73.7M 2s
    ##  70150K .......... .......... .......... .......... .......... 52%  113M 2s
    ##  70200K .......... .......... .......... .......... .......... 52% 7.00M 2s
    ##  70250K .......... .......... .......... .......... .......... 52%  121M 2s
    ##  70300K .......... .......... .......... .......... .......... 52%  130M 2s
    ##  70350K .......... .......... .......... .......... .......... 52% 11.9M 2s
    ##  70400K .......... .......... .......... .......... .......... 52%  112M 2s
    ##  70450K .......... .......... .......... .......... .......... 52%  139M 2s
    ##  70500K .......... .......... .......... .......... .......... 52% 14.7M 2s
    ##  70550K .......... .......... .......... .......... .......... 52% 86.3M 2s
    ##  70600K .......... .......... .......... .......... .......... 52%  136M 2s
    ##  70650K .......... .......... .......... .......... .......... 52% 13.3M 2s
    ##  70700K .......... .......... .......... .......... .......... 52%  130M 2s
    ##  70750K .......... .......... .......... .......... .......... 52% 13.9M 2s
    ##  70800K .......... .......... .......... .......... .......... 52%  114M 2s
    ##  70850K .......... .......... .......... .......... .......... 52%  117M 2s
    ##  70900K .......... .......... .......... .......... .......... 52% 13.2M 2s
    ##  70950K .......... .......... .......... .......... .......... 52%  104M 2s
    ##  71000K .......... .......... .......... .......... .......... 52% 14.3M 2s
    ##  71050K .......... .......... .......... .......... .......... 53% 84.3M 2s
    ##  71100K .......... .......... .......... .......... .......... 53%  132M 2s
    ##  71150K .......... .......... .......... .......... .......... 53% 14.4M 2s
    ##  71200K .......... .......... .......... .......... .......... 53%  114M 2s
    ##  71250K .......... .......... .......... .......... .......... 53% 15.5M 2s
    ##  71300K .......... .......... .......... .......... .......... 53% 84.7M 2s
    ##  71350K .......... .......... .......... .......... .......... 53%  115M 2s
    ##  71400K .......... .......... .......... .......... .......... 53% 12.4M 2s
    ##  71450K .......... .......... .......... .......... .......... 53%  123M 2s
    ##  71500K .......... .......... .......... .......... .......... 53%  132M 2s
    ##  71550K .......... .......... .......... .......... .......... 53% 12.3M 2s
    ##  71600K .......... .......... .......... .......... .......... 53%  107M 2s
    ##  71650K .......... .......... .......... .......... .......... 53% 12.2M 2s
    ##  71700K .......... .......... .......... .......... .......... 53%  111M 2s
    ##  71750K .......... .......... .......... .......... .......... 53%  124M 2s
    ##  71800K .......... .......... .......... .......... .......... 53% 11.8M 2s
    ##  71850K .......... .......... .......... .......... .......... 53%  130M 2s
    ##  71900K .......... .......... .......... .......... .......... 53% 12.3M 2s
    ##  71950K .......... .......... .......... .......... .......... 53% 82.5M 2s
    ##  72000K .......... .......... .......... .......... .......... 53%  128M 2s
    ##  72050K .......... .......... .......... .......... .......... 53% 14.9M 2s
    ##  72100K .......... .......... .......... .......... .......... 53%  107M 2s
    ##  72150K .......... .......... .......... .......... .......... 53% 14.6M 2s
    ##  72200K .......... .......... .......... .......... .......... 53%  102M 2s
    ##  72250K .......... .......... .......... .......... .......... 53%  107M 2s
    ##  72300K .......... .......... .......... .......... .......... 53% 16.4M 2s
    ##  72350K .......... .......... .......... .......... .......... 54% 99.4M 2s
    ##  72400K .......... .......... .......... .......... .......... 54% 90.2M 2s
    ##  72450K .......... .......... .......... .......... .......... 54% 14.9M 2s
    ##  72500K .......... .......... .......... .......... .......... 54% 98.3M 2s
    ##  72550K .......... .......... .......... .......... .......... 54%  101M 2s
    ##  72600K .......... .......... .......... .......... .......... 54% 78.2M 2s
    ##  72650K .......... .......... .......... .......... .......... 54% 18.3M 2s
    ##  72700K .......... .......... .......... .......... .......... 54%  106M 2s
    ##  72750K .......... .......... .......... .......... .......... 54% 46.0M 2s
    ##  72800K .......... .......... .......... .......... .......... 54% 23.0M 2s
    ##  72850K .......... .......... .......... .......... .......... 54% 28.8M 2s
    ##  72900K .......... .......... .......... .......... .......... 54% 23.2M 2s
    ##  72950K .......... .......... .......... .......... .......... 54%  130M 2s
    ##  73000K .......... .......... .......... .......... .......... 54% 30.0M 2s
    ##  73050K .......... .......... .......... .......... .......... 54% 22.5M 2s
    ##  73100K .......... .......... .......... .......... .......... 54% 36.2M 2s
    ##  73150K .......... .......... .......... .......... .......... 54% 19.9M 2s
    ##  73200K .......... .......... .......... .......... .......... 54%  140M 2s
    ##  73250K .......... .......... .......... .......... .......... 54% 42.7M 2s
    ##  73300K .......... .......... .......... .......... .......... 54% 10.1M 2s
    ##  73350K .......... .......... .......... .......... .......... 54%  124M 2s
    ##  73400K .......... .......... .......... .......... .......... 54%  164M 2s
    ##  73450K .......... .......... .......... .......... .......... 54%  169M 2s
    ##  73500K .......... .......... .......... .......... .......... 54% 13.8M 2s
    ##  73550K .......... .......... .......... .......... .......... 54%  105M 2s
    ##  73600K .......... .......... .......... .......... .......... 54% 19.0M 2s
    ##  73650K .......... .......... .......... .......... .......... 54% 49.6M 2s
    ##  73700K .......... .......... .......... .......... .......... 55% 57.7M 2s
    ##  73750K .......... .......... .......... .......... .......... 55% 19.9M 2s
    ##  73800K .......... .......... .......... .......... .......... 55% 95.2M 2s
    ##  73850K .......... .......... .......... .......... .......... 55%  105M 2s
    ##  73900K .......... .......... .......... .......... .......... 55% 12.9M 2s
    ##  73950K .......... .......... .......... .......... .......... 55% 87.3M 2s
    ##  74000K .......... .......... .......... .......... .......... 55% 21.8M 2s
    ##  74050K .......... .......... .......... .......... .......... 55% 34.5M 2s
    ##  74100K .......... .......... .......... .......... .......... 55%  147M 2s
    ##  74150K .......... .......... .......... .......... .......... 55% 12.5M 2s
    ##  74200K .......... .......... .......... .......... .......... 55%  106M 2s
    ##  74250K .......... .......... .......... .......... .......... 55% 14.5M 2s
    ##  74300K .......... .......... .......... .......... .......... 55%  125M 2s
    ##  74350K .......... .......... .......... .......... .......... 55% 13.2M 2s
    ##  74400K .......... .......... .......... .......... .......... 55%  124M 2s
    ##  74450K .......... .......... .......... .......... .......... 55%  137M 2s
    ##  74500K .......... .......... .......... .......... .......... 55% 13.4M 2s
    ##  74550K .......... .......... .......... .......... .......... 55%  131M 2s
    ##  74600K .......... .......... .......... .......... .......... 55%  167M 2s
    ##  74650K .......... .......... .......... .......... .......... 55% 14.8M 2s
    ##  74700K .......... .......... .......... .......... .......... 55%  123M 2s
    ##  74750K .......... .......... .......... .......... .......... 55% 14.0M 2s
    ##  74800K .......... .......... .......... .......... .......... 55%  128M 2s
    ##  74850K .......... .......... .......... .......... .......... 55%  168M 2s
    ##  74900K .......... .......... .......... .......... .......... 55% 15.3M 2s
    ##  74950K .......... .......... .......... .......... .......... 55%  124M 2s
    ##  75000K .......... .......... .......... .......... .......... 55% 14.2M 2s
    ##  75050K .......... .......... .......... .......... .......... 56%  135M 2s
    ##  75100K .......... .......... .......... .......... .......... 56%  158M 2s
    ##  75150K .......... .......... .......... .......... .......... 56% 12.9M 2s
    ##  75200K .......... .......... .......... .......... .......... 56%  130M 2s
    ##  75250K .......... .......... .......... .......... .......... 56% 12.7M 2s
    ##  75300K .......... .......... .......... .......... .......... 56%  122M 2s
    ##  75350K .......... .......... .......... .......... .......... 56%  139M 2s
    ##  75400K .......... .......... .......... .......... .......... 56% 12.4M 2s
    ##  75450K .......... .......... .......... .......... .......... 56%  110M 2s
    ##  75500K .......... .......... .......... .......... .......... 56%  157M 2s
    ##  75550K .......... .......... .......... .......... .......... 56% 14.0M 2s
    ##  75600K .......... .......... .......... .......... .......... 56%  116M 2s
    ##  75650K .......... .......... .......... .......... .......... 56% 15.2M 2s
    ##  75700K .......... .......... .......... .......... .......... 56% 87.5M 2s
    ##  75750K .......... .......... .......... .......... .......... 56%  133M 2s
    ##  75800K .......... .......... .......... .......... .......... 56% 14.6M 2s
    ##  75850K .......... .......... .......... .......... .......... 56% 60.6M 2s
    ##  75900K .......... .......... .......... .......... .......... 56% 18.3M 2s
    ##  75950K .......... .......... .......... .......... .......... 56% 40.6M 2s
    ##  76000K .......... .......... .......... .......... .......... 56%  133M 2s
    ##  76050K .......... .......... .......... .......... .......... 56% 17.0M 2s
    ##  76100K .......... .......... .......... .......... .......... 56% 96.5M 2s
    ##  76150K .......... .......... .......... .......... .......... 56% 15.4M 2s
    ##  76200K .......... .......... .......... .......... .......... 56% 55.9M 2s
    ##  76250K .......... .......... .......... .......... .......... 56%  126M 2s
    ##  76300K .......... .......... .......... .......... .......... 56% 15.8M 2s
    ##  76350K .......... .......... .......... .......... .......... 56% 83.3M 2s
    ##  76400K .......... .......... .......... .......... .......... 57% 14.3M 2s
    ##  76450K .......... .......... .......... .......... .......... 57% 87.8M 2s
    ##  76500K .......... .......... .......... .......... .......... 57%  127M 2s
    ##  76550K .......... .......... .......... .......... .......... 57% 17.7M 2s
    ##  76600K .......... .......... .......... .......... .......... 57% 48.6M 2s
    ##  76650K .......... .......... .......... .......... .......... 57% 19.0M 2s
    ##  76700K .......... .......... .......... .......... .......... 57%  145M 2s
    ##  76750K .......... .......... .......... .......... .......... 57% 59.2M 2s
    ##  76800K .......... .......... .......... .......... .......... 57% 19.3M 2s
    ##  76850K .......... .......... .......... .......... .......... 57% 43.0M 2s
    ##  76900K .......... .......... .......... .......... .......... 57% 20.0M 2s
    ##  76950K .......... .......... .......... .......... .......... 57% 97.3M 2s
    ##  77000K .......... .......... .......... .......... .......... 57% 46.8M 2s
    ##  77050K .......... .......... .......... .......... .......... 57% 18.4M 2s
    ##  77100K .......... .......... .......... .......... .......... 57% 41.9M 2s
    ##  77150K .......... .......... .......... .......... .......... 57% 28.3M 2s
    ##  77200K .......... .......... .......... .......... .......... 57%  113M 2s
    ##  77250K .......... .......... .......... .......... .......... 57% 25.8M 2s
    ##  77300K .......... .......... .......... .......... .......... 57% 27.0M 2s
    ##  77350K .......... .......... .......... .......... .......... 57% 33.4M 2s
    ##  77400K .......... .......... .......... .......... .......... 57% 88.5M 2s
    ##  77450K .......... .......... .......... .......... .......... 57% 20.6M 2s
    ##  77500K .......... .......... .......... .......... .......... 57% 35.0M 2s
    ##  77550K .......... .......... .......... .......... .......... 57% 18.2M 2s
    ##  77600K .......... .......... .......... .......... .......... 57%  118M 2s
    ##  77650K .......... .......... .......... .......... .......... 57% 53.5M 2s
    ##  77700K .......... .......... .......... .......... .......... 57% 16.7M 2s
    ##  77750K .......... .......... .......... .......... .......... 58% 80.7M 2s
    ##  77800K .......... .......... .......... .......... .......... 58% 13.5M 2s
    ##  77850K .......... .......... .......... .......... .......... 58% 89.2M 2s
    ##  77900K .......... .......... .......... .......... .......... 58%  145M 2s
    ##  77950K .......... .......... .......... .......... .......... 58% 14.4M 2s
    ##  78000K .......... .......... .......... .......... .......... 58% 44.1M 2s
    ##  78050K .......... .......... .......... .......... .......... 58% 75.7M 2s
    ##  78100K .......... .......... .......... .......... .......... 58% 85.0M 2s
    ##  78150K .......... .......... .......... .......... .......... 58% 6.83M 2s
    ##  78200K .......... .......... .......... .......... .......... 58%  147M 2s
    ##  78250K .......... .......... .......... .......... .......... 58%  166M 2s
    ##  78300K .......... .......... .......... .......... .......... 58%  152M 2s
    ##  78350K .......... .......... .......... .......... .......... 58% 14.2M 2s
    ##  78400K .......... .......... .......... .......... .......... 58% 80.3M 2s
    ##  78450K .......... .......... .......... .......... .......... 58% 79.9M 2s
    ##  78500K .......... .......... .......... .......... .......... 58% 17.5M 2s
    ##  78550K .......... .......... .......... .......... .......... 58% 78.8M 2s
    ##  78600K .......... .......... .......... .......... .......... 58%  156M 2s
    ##  78650K .......... .......... .......... .......... .......... 58% 52.7M 2s
    ##  78700K .......... .......... .......... .......... .......... 58% 24.1M 2s
    ##  78750K .......... .......... .......... .......... .......... 58% 86.8M 2s
    ##  78800K .......... .......... .......... .......... .......... 58% 38.0M 2s
    ##  78850K .......... .......... .......... .......... .......... 58%  123M 2s
    ##  78900K .......... .......... .......... .......... .......... 58% 29.2M 2s
    ##  78950K .......... .......... .......... .......... .......... 58% 23.0M 2s
    ##  79000K .......... .......... .......... .......... .......... 58%  143M 2s
    ##  79050K .......... .......... .......... .......... .......... 59%  147M 2s
    ##  79100K .......... .......... .......... .......... .......... 59% 48.0M 2s
    ##  79150K .......... .......... .......... .......... .......... 59% 22.1M 2s
    ##  79200K .......... .......... .......... .......... .......... 59%  126M 2s
    ##  79250K .......... .......... .......... .......... .......... 59%  128M 2s
    ##  79300K .......... .......... .......... .......... .......... 59%  118M 2s
    ##  79350K .......... .......... .......... .......... .......... 59% 18.0M 2s
    ##  79400K .......... .......... .......... .......... .......... 59%  149M 2s
    ##  79450K .......... .......... .......... .......... .......... 59% 60.6M 2s
    ##  79500K .......... .......... .......... .......... .......... 59%  145M 2s
    ##  79550K .......... .......... .......... .......... .......... 59% 18.9M 2s
    ##  79600K .......... .......... .......... .......... .......... 59% 97.6M 2s
    ##  79650K .......... .......... .......... .......... .......... 59% 57.4M 2s
    ##  79700K .......... .......... .......... .......... .......... 59% 24.5M 2s
    ##  79750K .......... .......... .......... .......... .......... 59%  110M 2s
    ##  79800K .......... .......... .......... .......... .......... 59% 27.7M 2s
    ##  79850K .......... .......... .......... .......... .......... 59%  127M 2s
    ##  79900K .......... .......... .......... .......... .......... 59% 41.8M 2s
    ##  79950K .......... .......... .......... .......... .......... 59% 27.9M 2s
    ##  80000K .......... .......... .......... .......... .......... 59%  105M 2s
    ##  80050K .......... .......... .......... .......... .......... 59% 98.8M 2s
    ##  80100K .......... .......... .......... .......... .......... 59% 71.2M 2s
    ##  80150K .......... .......... .......... .......... .......... 59% 31.1M 2s
    ##  80200K .......... .......... .......... .......... .......... 59% 35.2M 2s
    ##  80250K .......... .......... .......... .......... .......... 59%  126M 2s
    ##  80300K .......... .......... .......... .......... .......... 59% 44.3M 2s
    ##  80350K .......... .......... .......... .......... .......... 59% 21.3M 2s
    ##  80400K .......... .......... .......... .......... .......... 60% 34.5M 2s
    ##  80450K .......... .......... .......... .......... .......... 60% 96.7M 2s
    ##  80500K .......... .......... .......... .......... .......... 60%  128M 2s
    ##  80550K .......... .......... .......... .......... .......... 60% 25.9M 2s
    ##  80600K .......... .......... .......... .......... .......... 60% 27.6M 2s
    ##  80650K .......... .......... .......... .......... .......... 60%  114M 2s
    ##  80700K .......... .......... .......... .......... .......... 60%  151M 2s
    ##  80750K .......... .......... .......... .......... .......... 60% 17.5M 2s
    ##  80800K .......... .......... .......... .......... .......... 60% 38.4M 2s
    ##  80850K .......... .......... .......... .......... .......... 60%  117M 2s
    ##  80900K .......... .......... .......... .......... .......... 60% 23.8M 2s
    ##  80950K .......... .......... .......... .......... .......... 60%  130M 2s
    ##  81000K .......... .......... .......... .......... .......... 60% 19.6M 2s
    ##  81050K .......... .......... .......... .......... .......... 60%  118M 2s
    ##  81100K .......... .......... .......... .......... .......... 60%  119M 2s
    ##  81150K .......... .......... .......... .......... .......... 60% 19.9M 2s
    ##  81200K .......... .......... .......... .......... .......... 60% 72.3M 2s
    ##  81250K .......... .......... .......... .......... .......... 60% 80.2M 2s
    ##  81300K .......... .......... .......... .......... .......... 60% 56.8M 2s
    ##  81350K .......... .......... .......... .......... .......... 60% 45.2M 2s
    ##  81400K .......... .......... .......... .......... .......... 60% 33.5M 2s
    ##  81450K .......... .......... .......... .......... .......... 60% 95.3M 2s
    ##  81500K .......... .......... .......... .......... .......... 60% 71.7M 2s
    ##  81550K .......... .......... .......... .......... .......... 60% 13.2M 2s
    ##  81600K .......... .......... .......... .......... .......... 60%  117M 2s
    ##  81650K .......... .......... .......... .......... .......... 60%  115M 2s
    ##  81700K .......... .......... .......... .......... .......... 60%  162M 2s
    ##  81750K .......... .......... .......... .......... .......... 61% 12.1M 2s
    ##  81800K .......... .......... .......... .......... .......... 61% 98.0M 2s
    ##  81850K .......... .......... .......... .......... .......... 61%  153M 2s
    ##  81900K .......... .......... .......... .......... .......... 61%  128M 2s
    ##  81950K .......... .......... .......... .......... .......... 61% 15.6M 2s
    ##  82000K .......... .......... .......... .......... .......... 61% 98.4M 2s
    ##  82050K .......... .......... .......... .......... .......... 61%  127M 2s
    ##  82100K .......... .......... .......... .......... .......... 61%  121M 2s
    ##  82150K .......... .......... .......... .......... .......... 61% 6.61M 2s
    ##  82200K .......... .......... .......... .......... .......... 61%  109M 2s
    ##  82250K .......... .......... .......... .......... .......... 61%  156M 2s
    ##  82300K .......... .......... .......... .......... .......... 61%  163M 2s
    ##  82350K .......... .......... .......... .......... .......... 61% 11.2M 2s
    ##  82400K .......... .......... .......... .......... .......... 61%  127M 2s
    ##  82450K .......... .......... .......... .......... .......... 61%  161M 2s
    ##  82500K .......... .......... .......... .......... .......... 61% 16.8M 2s
    ##  82550K .......... .......... .......... .......... .......... 61% 91.9M 2s
    ##  82600K .......... .......... .......... .......... .......... 61% 43.1M 2s
    ##  82650K .......... .......... .......... .......... .......... 61%  108M 2s
    ##  82700K .......... .......... .......... .......... .......... 61% 24.0M 2s
    ##  82750K .......... .......... .......... .......... .......... 61% 26.4M 2s
    ##  82800K .......... .......... .......... .......... .......... 61%  117M 2s
    ##  82850K .......... .......... .......... .......... .......... 61%  153M 2s
    ##  82900K .......... .......... .......... .......... .......... 61% 19.2M 2s
    ##  82950K .......... .......... .......... .......... .......... 61% 42.0M 2s
    ##  83000K .......... .......... .......... .......... .......... 61%  126M 2s
    ##  83050K .......... .......... .......... .......... .......... 61%  126M 2s
    ##  83100K .......... .......... .......... .......... .......... 62% 24.7M 2s
    ##  83150K .......... .......... .......... .......... .......... 62% 29.5M 2s
    ##  83200K .......... .......... .......... .......... .......... 62% 98.9M 2s
    ##  83250K .......... .......... .......... .......... .......... 62% 28.7M 2s
    ##  83300K .......... .......... .......... .......... .......... 62% 88.4M 2s
    ##  83350K .......... .......... .......... .......... .......... 62% 23.5M 2s
    ##  83400K .......... .......... .......... .......... .......... 62% 83.5M 2s
    ##  83450K .......... .......... .......... .......... .......... 62% 72.6M 2s
    ##  83500K .......... .......... .......... .......... .......... 62% 95.0M 2s
    ##  83550K .......... .......... .......... .......... .......... 62% 21.5M 2s
    ##  83600K .......... .......... .......... .......... .......... 62% 37.3M 2s
    ##  83650K .......... .......... .......... .......... .......... 62% 74.6M 2s
    ##  83700K .......... .......... .......... .......... .......... 62% 45.1M 2s
    ##  83750K .......... .......... .......... .......... .......... 62% 76.3M 2s
    ##  83800K .......... .......... .......... .......... .......... 62% 28.6M 2s
    ##  83850K .......... .......... .......... .......... .......... 62% 67.2M 2s
    ##  83900K .......... .......... .......... .......... .......... 62% 73.4M 2s
    ##  83950K .......... .......... .......... .......... .......... 62% 27.5M 2s
    ##  84000K .......... .......... .......... .......... .......... 62% 83.3M 2s
    ##  84050K .......... .......... .......... .......... .......... 62% 83.1M 2s
    ##  84100K .......... .......... .......... .......... .......... 62% 24.8M 2s
    ##  84150K .......... .......... .......... .......... .......... 62% 19.4M 2s
    ##  84200K .......... .......... .......... .......... .......... 62% 89.5M 2s
    ##  84250K .......... .......... .......... .......... .......... 62%  111M 2s
    ##  84300K .......... .......... .......... .......... .......... 62% 15.6M 2s
    ##  84350K .......... .......... .......... .......... .......... 62% 60.6M 2s
    ##  84400K .......... .......... .......... .......... .......... 62% 74.9M 2s
    ##  84450K .......... .......... .......... .......... .......... 63%  113M 2s
    ##  84500K .......... .......... .......... .......... .......... 63% 28.0M 2s
    ##  84550K .......... .......... .......... .......... .......... 63% 27.9M 2s
    ##  84600K .......... .......... .......... .......... .......... 63% 96.6M 2s
    ##  84650K .......... .......... .......... .......... .......... 63% 58.8M 2s
    ##  84700K .......... .......... .......... .......... .......... 63% 49.8M 2s
    ##  84750K .......... .......... .......... .......... .......... 63% 24.0M 2s
    ##  84800K .......... .......... .......... .......... .......... 63% 97.2M 2s
    ##  84850K .......... .......... .......... .......... .......... 63% 86.2M 2s
    ##  84900K .......... .......... .......... .......... .......... 63% 10.3M 2s
    ##  84950K .......... .......... .......... .......... .......... 63% 85.1M 2s
    ##  85000K .......... .......... .......... .......... .......... 63% 91.5M 2s
    ##  85050K .......... .......... .......... .......... .......... 63%  127M 2s
    ##  85100K .......... .......... .......... .......... .......... 63% 15.3M 2s
    ##  85150K .......... .......... .......... .......... .......... 63% 84.4M 2s
    ##  85200K .......... .......... .......... .......... .......... 63% 99.7M 2s
    ##  85250K .......... .......... .......... .......... .......... 63% 96.9M 2s
    ##  85300K .......... .......... .......... .......... .......... 63% 18.4M 2s
    ##  85350K .......... .......... .......... .......... .......... 63% 89.7M 2s
    ##  85400K .......... .......... .......... .......... .......... 63% 94.1M 2s
    ##  85450K .......... .......... .......... .......... .......... 63% 83.0M 2s
    ##  85500K .......... .......... .......... .......... .......... 63% 12.6M 2s
    ##  85550K .......... .......... .......... .......... .......... 63% 71.0M 2s
    ##  85600K .......... .......... .......... .......... .......... 63%  123M 2s
    ##  85650K .......... .......... .......... .......... .......... 63%  119M 2s
    ##  85700K .......... .......... .......... .......... .......... 63% 13.6M 2s
    ##  85750K .......... .......... .......... .......... .......... 63% 72.9M 2s
    ##  85800K .......... .......... .......... .......... .......... 64% 93.8M 2s
    ##  85850K .......... .......... .......... .......... .......... 64%  123M 2s
    ##  85900K .......... .......... .......... .......... .......... 64% 22.6M 2s
    ##  85950K .......... .......... .......... .......... .......... 64% 76.1M 2s
    ##  86000K .......... .......... .......... .......... .......... 64% 54.8M 2s
    ##  86050K .......... .......... .......... .......... .......... 64% 23.1M 2s
    ##  86100K .......... .......... .......... .......... .......... 64% 90.0M 2s
    ##  86150K .......... .......... .......... .......... .......... 64% 99.8M 2s
    ##  86200K .......... .......... .......... .......... .......... 64% 49.3M 2s
    ##  86250K .......... .......... .......... .......... .......... 64% 22.6M 2s
    ##  86300K .......... .......... .......... .......... .......... 64% 85.3M 2s
    ##  86350K .......... .......... .......... .......... .......... 64% 92.1M 2s
    ##  86400K .......... .......... .......... .......... .......... 64% 84.1M 2s
    ##  86450K .......... .......... .......... .......... .......... 64% 17.4M 2s
    ##  86500K .......... .......... .......... .......... .......... 64% 97.1M 2s
    ##  86550K .......... .......... .......... .......... .......... 64% 84.8M 2s
    ##  86600K .......... .......... .......... .......... .......... 64%  122M 2s
    ##  86650K .......... .......... .......... .......... .......... 64% 24.1M 2s
    ##  86700K .......... .......... .......... .......... .......... 64% 51.5M 2s
    ##  86750K .......... .......... .......... .......... .......... 64% 66.9M 2s
    ##  86800K .......... .......... .......... .......... .......... 64% 97.7M 2s
    ##  86850K .......... .......... .......... .......... .......... 64% 49.6M 2s
    ##  86900K .......... .......... .......... .......... .......... 64% 26.7M 2s
    ##  86950K .......... .......... .......... .......... .......... 64% 76.5M 2s
    ##  87000K .......... .......... .......... .......... .......... 64% 57.6M 2s
    ##  87050K .......... .......... .......... .......... .......... 64% 94.9M 2s
    ##  87100K .......... .......... .......... .......... .......... 65% 21.0M 1s
    ##  87150K .......... .......... .......... .......... .......... 65% 97.8M 1s
    ##  87200K .......... .......... .......... .......... .......... 65%  111M 1s
    ##  87250K .......... .......... .......... .......... .......... 65% 18.4M 1s
    ##  87300K .......... .......... .......... .......... .......... 65%  114M 1s
    ##  87350K .......... .......... .......... .......... .......... 65% 59.2M 1s
    ##  87400K .......... .......... .......... .......... .......... 65% 88.2M 1s
    ##  87450K .......... .......... .......... .......... .......... 65% 24.7M 1s
    ##  87500K .......... .......... .......... .......... .......... 65% 97.0M 1s
    ##  87550K .......... .......... .......... .......... .......... 65% 33.6M 1s
    ##  87600K .......... .......... .......... .......... .......... 65% 21.6M 1s
    ##  87650K .......... .......... .......... .......... .......... 65%  127M 1s
    ##  87700K .......... .......... .......... .......... .......... 65% 61.8M 1s
    ##  87750K .......... .......... .......... .......... .......... 65%  107M 1s
    ##  87800K .......... .......... .......... .......... .......... 65% 35.9M 1s
    ##  87850K .......... .......... .......... .......... .......... 65% 20.1M 1s
    ##  87900K .......... .......... .......... .......... .......... 65%  116M 1s
    ##  87950K .......... .......... .......... .......... .......... 65% 54.1M 1s
    ##  88000K .......... .......... .......... .......... .......... 65% 95.3M 1s
    ##  88050K .......... .......... .......... .......... .......... 65% 28.8M 1s
    ##  88100K .......... .......... .......... .......... .......... 65% 54.4M 1s
    ##  88150K .......... .......... .......... .......... .......... 65% 40.0M 1s
    ##  88200K .......... .......... .......... .......... .......... 65%  141M 1s
    ##  88250K .......... .......... .......... .......... .......... 65% 55.2M 1s
    ##  88300K .......... .......... .......... .......... .......... 65% 27.4M 1s
    ##  88350K .......... .......... .......... .......... .......... 65% 62.1M 1s
    ##  88400K .......... .......... .......... .......... .......... 65% 40.8M 1s
    ##  88450K .......... .......... .......... .......... .......... 66%  124M 1s
    ##  88500K .......... .......... .......... .......... .......... 66% 39.5M 1s
    ##  88550K .......... .......... .......... .......... .......... 66% 31.5M 1s
    ##  88600K .......... .......... .......... .......... .......... 66% 53.4M 1s
    ##  88650K .......... .......... .......... .......... .......... 66%  143M 1s
    ##  88700K .......... .......... .......... .......... .......... 66% 40.2M 1s
    ##  88750K .......... .......... .......... .......... .......... 66% 22.0M 1s
    ##  88800K .......... .......... .......... .......... .......... 66%  130M 1s
    ##  88850K .......... .......... .......... .......... .......... 66% 50.7M 1s
    ##  88900K .......... .......... .......... .......... .......... 66%  132M 1s
    ##  88950K .......... .......... .......... .......... .......... 66% 18.5M 1s
    ##  89000K .......... .......... .......... .......... .......... 66%  120M 1s
    ##  89050K .......... .......... .......... .......... .......... 66%  126M 1s
    ##  89100K .......... .......... .......... .......... .......... 66% 16.2M 1s
    ##  89150K .......... .......... .......... .......... .......... 66% 75.0M 1s
    ##  89200K .......... .......... .......... .......... .......... 66%  129M 1s
    ##  89250K .......... .......... .......... .......... .......... 66% 70.7M 1s
    ##  89300K .......... .......... .......... .......... .......... 66% 17.1M 1s
    ##  89350K .......... .......... .......... .......... .......... 66%  110M 1s
    ##  89400K .......... .......... .......... .......... .......... 66%  142M 1s
    ##  89450K .......... .......... .......... .......... .......... 66%  112M 1s
    ##  89500K .......... .......... .......... .......... .......... 66% 15.9M 1s
    ##  89550K .......... .......... .......... .......... .......... 66% 88.7M 1s
    ##  89600K .......... .......... .......... .......... .......... 66%  125M 1s
    ##  89650K .......... .......... .......... .......... .......... 66%  111M 1s
    ##  89700K .......... .......... .......... .......... .......... 66% 15.8M 1s
    ##  89750K .......... .......... .......... .......... .......... 66%  106M 1s
    ##  89800K .......... .......... .......... .......... .......... 67%  102M 1s
    ##  89850K .......... .......... .......... .......... .......... 67%  110M 1s
    ##  89900K .......... .......... .......... .......... .......... 67% 16.4M 1s
    ##  89950K .......... .......... .......... .......... .......... 67%  102M 1s
    ##  90000K .......... .......... .......... .......... .......... 67% 92.0M 1s
    ##  90050K .......... .......... .......... .......... .......... 67% 42.8M 1s
    ##  90100K .......... .......... .......... .......... .......... 67% 48.8M 1s
    ##  90150K .......... .......... .......... .......... .......... 67% 89.5M 1s
    ##  90200K .......... .......... .......... .......... .......... 67% 28.1M 1s
    ##  90250K .......... .......... .......... .......... .......... 67% 65.5M 1s
    ##  90300K .......... .......... .......... .......... .......... 67%  127M 1s
    ##  90350K .......... .......... .......... .......... .......... 67% 78.6M 1s
    ##  90400K .......... .......... .......... .......... .......... 67% 20.9M 1s
    ##  90450K .......... .......... .......... .......... .......... 67%  104M 1s
    ##  90500K .......... .......... .......... .......... .......... 67% 48.4M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 95.7M 1s
    ##  90600K .......... .......... .......... .......... .......... 67% 47.2M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 32.1M 1s
    ##  90700K .......... .......... .......... .......... .......... 67% 33.7M 1s
    ##  90750K .......... .......... .......... .......... .......... 67% 78.2M 1s
    ##  90800K .......... .......... .......... .......... .......... 67%  141M 1s
    ##  90850K .......... .......... .......... .......... .......... 67% 16.3M 1s
    ##  90900K .......... .......... .......... .......... .......... 67% 85.4M 1s
    ##  90950K .......... .......... .......... .......... .......... 67% 97.1M 1s
    ##  91000K .......... .......... .......... .......... .......... 67%  116M 1s
    ##  91050K .......... .......... .......... .......... .......... 67% 22.2M 1s
    ##  91100K .......... .......... .......... .......... .......... 67% 94.4M 1s
    ##  91150K .......... .......... .......... .......... .......... 68% 40.1M 1s
    ##  91200K .......... .......... .......... .......... .......... 68% 49.7M 1s
    ##  91250K .......... .......... .......... .......... .......... 68% 72.2M 1s
    ##  91300K .......... .......... .......... .......... .......... 68% 71.2M 1s
    ##  91350K .......... .......... .......... .......... .......... 68% 42.0M 1s
    ##  91400K .......... .......... .......... .......... .......... 68% 61.8M 1s
    ##  91450K .......... .......... .......... .......... .......... 68% 30.1M 1s
    ##  91500K .......... .......... .......... .......... .......... 68%  113M 1s
    ##  91550K .......... .......... .......... .......... .......... 68% 31.8M 1s
    ##  91600K .......... .......... .......... .......... .......... 68%  105M 1s
    ##  91650K .......... .......... .......... .......... .......... 68% 32.3M 1s
    ##  91700K .......... .......... .......... .......... .......... 68% 26.7M 1s
    ##  91750K .......... .......... .......... .......... .......... 68% 63.6M 1s
    ##  91800K .......... .......... .......... .......... .......... 68%  138M 1s
    ##  91850K .......... .......... .......... .......... .......... 68%  107M 1s
    ##  91900K .......... .......... .......... .......... .......... 68% 20.8M 1s
    ##  91950K .......... .......... .......... .......... .......... 68% 71.9M 1s
    ##  92000K .......... .......... .......... .......... .......... 68%  117M 1s
    ##  92050K .......... .......... .......... .......... .......... 68% 59.3M 1s
    ##  92100K .......... .......... .......... .......... .......... 68% 20.2M 1s
    ##  92150K .......... .......... .......... .......... .......... 68% 95.8M 1s
    ##  92200K .......... .......... .......... .......... .......... 68%  128M 1s
    ##  92250K .......... .......... .......... .......... .......... 68% 94.1M 1s
    ##  92300K .......... .......... .......... .......... .......... 68% 18.1M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 78.1M 1s
    ##  92400K .......... .......... .......... .......... .......... 68%  125M 1s
    ##  92450K .......... .......... .......... .......... .......... 68%  141M 1s
    ##  92500K .......... .......... .......... .......... .......... 69% 16.7M 1s
    ##  92550K .......... .......... .......... .......... .......... 69% 83.5M 1s
    ##  92600K .......... .......... .......... .......... .......... 69%  130M 1s
    ##  92650K .......... .......... .......... .......... .......... 69%  137M 1s
    ##  92700K .......... .......... .......... .......... .......... 69% 20.9M 1s
    ##  92750K .......... .......... .......... .......... .......... 69% 81.9M 1s
    ##  92800K .......... .......... .......... .......... .......... 69%  117M 1s
    ##  92850K .......... .......... .......... .......... .......... 69% 25.9M 1s
    ##  92900K .......... .......... .......... .......... .......... 69%  100M 1s
    ##  92950K .......... .......... .......... .......... .......... 69% 96.1M 1s
    ##  93000K .......... .......... .......... .......... .......... 69% 89.2M 1s
    ##  93050K .......... .......... .......... .......... .......... 69% 22.2M 1s
    ##  93100K .......... .......... .......... .......... .......... 69% 83.2M 1s
    ##  93150K .......... .......... .......... .......... .......... 69%  111M 1s
    ##  93200K .......... .......... .......... .......... .......... 69%  110M 1s
    ##  93250K .......... .......... .......... .......... .......... 69% 22.6M 1s
    ##  93300K .......... .......... .......... .......... .......... 69% 86.2M 1s
    ##  93350K .......... .......... .......... .......... .......... 69%  112M 1s
    ##  93400K .......... .......... .......... .......... .......... 69%  121M 1s
    ##  93450K .......... .......... .......... .......... .......... 69% 18.1M 1s
    ##  93500K .......... .......... .......... .......... .......... 69% 91.2M 1s
    ##  93550K .......... .......... .......... .......... .......... 69% 52.5M 1s
    ##  93600K .......... .......... .......... .......... .......... 69% 20.9M 1s
    ##  93650K .......... .......... .......... .......... .......... 69% 91.5M 1s
    ##  93700K .......... .......... .......... .......... .......... 69%  124M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 79.2M 1s
    ##  93800K .......... .......... .......... .......... .......... 70% 17.3M 1s
    ##  93850K .......... .......... .......... .......... .......... 70%  109M 1s
    ##  93900K .......... .......... .......... .......... .......... 70%  104M 1s
    ##  93950K .......... .......... .......... .......... .......... 70% 68.4M 1s
    ##  94000K .......... .......... .......... .......... .......... 70% 18.4M 1s
    ##  94050K .......... .......... .......... .......... .......... 70%  103M 1s
    ##  94100K .......... .......... .......... .......... .......... 70%  133M 1s
    ##  94150K .......... .......... .......... .......... .......... 70% 92.2M 1s
    ##  94200K .......... .......... .......... .......... .......... 70% 14.2M 1s
    ##  94250K .......... .......... .......... .......... .......... 70% 59.5M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 89.0M 1s
    ##  94350K .......... .......... .......... .......... .......... 70% 41.5M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 25.8M 1s
    ##  94450K .......... .......... .......... .......... .......... 70% 94.8M 1s
    ##  94500K .......... .......... .......... .......... .......... 70% 53.1M 1s
    ##  94550K .......... .......... .......... .......... .......... 70% 79.1M 1s
    ##  94600K .......... .......... .......... .......... .......... 70% 25.2M 1s
    ##  94650K .......... .......... .......... .......... .......... 70% 26.2M 1s
    ##  94700K .......... .......... .......... .......... .......... 70%  121M 1s
    ##  94750K .......... .......... .......... .......... .......... 70% 38.2M 1s
    ##  94800K .......... .......... .......... .......... .......... 70% 86.1M 1s
    ##  94850K .......... .......... .......... .......... .......... 70% 20.0M 1s
    ##  94900K .......... .......... .......... .......... .......... 70%  112M 1s
    ##  94950K .......... .......... .......... .......... .......... 70%  116M 1s
    ##  95000K .......... .......... .......... .......... .......... 70%  135M 1s
    ##  95050K .......... .......... .......... .......... .......... 70% 12.5M 1s
    ##  95100K .......... .......... .......... .......... .......... 70% 99.4M 1s
    ##  95150K .......... .......... .......... .......... .......... 71%  109M 1s
    ##  95200K .......... .......... .......... .......... .......... 71% 14.3M 1s
    ##  95250K .......... .......... .......... .......... .......... 71% 88.3M 1s
    ##  95300K .......... .......... .......... .......... .......... 71% 89.1M 1s
    ##  95350K .......... .......... .......... .......... .......... 71%  115M 1s
    ##  95400K .......... .......... .......... .......... .......... 71% 19.2M 1s
    ##  95450K .......... .......... .......... .......... .......... 71% 97.0M 1s
    ##  95500K .......... .......... .......... .......... .......... 71%  106M 1s
    ##  95550K .......... .......... .......... .......... .......... 71% 81.9M 1s
    ##  95600K .......... .......... .......... .......... .......... 71% 17.1M 1s
    ##  95650K .......... .......... .......... .......... .......... 71%  100M 1s
    ##  95700K .......... .......... .......... .......... .......... 71% 56.1M 1s
    ##  95750K .......... .......... .......... .......... .......... 71% 82.1M 1s
    ##  95800K .......... .......... .......... .......... .......... 71% 22.8M 1s
    ##  95850K .......... .......... .......... .......... .......... 71%  110M 1s
    ##  95900K .......... .......... .......... .......... .......... 71% 65.7M 1s
    ##  95950K .......... .......... .......... .......... .......... 71% 30.9M 1s
    ##  96000K .......... .......... .......... .......... .......... 71% 91.2M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 25.8M 1s
    ##  96100K .......... .......... .......... .......... .......... 71%  112M 1s
    ##  96150K .......... .......... .......... .......... .......... 71% 40.6M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 66.1M 1s
    ##  96250K .......... .......... .......... .......... .......... 71%  109M 1s
    ##  96300K .......... .......... .......... .......... .......... 71% 80.3M 1s
    ##  96350K .......... .......... .......... .......... .......... 71% 20.9M 1s
    ##  96400K .......... .......... .......... .......... .......... 71% 74.7M 1s
    ##  96450K .......... .......... .......... .......... .......... 71% 89.3M 1s
    ##  96500K .......... .......... .......... .......... .......... 72% 94.2M 1s
    ##  96550K .......... .......... .......... .......... .......... 72% 39.4M 1s
    ##  96600K .......... .......... .......... .......... .......... 72% 50.2M 1s
    ##  96650K .......... .......... .......... .......... .......... 72% 92.2M 1s
    ##  96700K .......... .......... .......... .......... .......... 72% 9.19M 1s
    ##  96750K .......... .......... .......... .......... .......... 72% 62.3M 1s
    ##  96800K .......... .......... .......... .......... .......... 72% 87.6M 1s
    ##  96850K .......... .......... .......... .......... .......... 72%  104M 1s
    ##  96900K .......... .......... .......... .......... .......... 72%  119M 1s
    ##  96950K .......... .......... .......... .......... .......... 72% 14.9M 1s
    ##  97000K .......... .......... .......... .......... .......... 72% 96.4M 1s
    ##  97050K .......... .......... .......... .......... .......... 72%  105M 1s
    ##  97100K .......... .......... .......... .......... .......... 72%  108M 1s
    ##  97150K .......... .......... .......... .......... .......... 72%  102M 1s
    ##  97200K .......... .......... .......... .......... .......... 72% 16.4M 1s
    ##  97250K .......... .......... .......... .......... .......... 72% 38.8M 1s
    ##  97300K .......... .......... .......... .......... .......... 72% 73.5M 1s
    ##  97350K .......... .......... .......... .......... .......... 72% 88.8M 1s
    ##  97400K .......... .......... .......... .......... .......... 72%  111M 1s
    ##  97450K .......... .......... .......... .......... .......... 72% 45.7M 1s
    ##  97500K .......... .......... .......... .......... .......... 72% 19.8M 1s
    ##  97550K .......... .......... .......... .......... .......... 72% 91.2M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 98.3M 1s
    ##  97650K .......... .......... .......... .......... .......... 72%  101M 1s
    ##  97700K .......... .......... .......... .......... .......... 72% 55.4M 1s
    ##  97750K .......... .......... .......... .......... .......... 72% 47.8M 1s
    ##  97800K .......... .......... .......... .......... .......... 72% 26.0M 1s
    ##  97850K .......... .......... .......... .......... .......... 73% 78.4M 1s
    ##  97900K .......... .......... .......... .......... .......... 73%  129M 1s
    ##  97950K .......... .......... .......... .......... .......... 73% 64.3M 1s
    ##  98000K .......... .......... .......... .......... .......... 73%  101M 1s
    ##  98050K .......... .......... .......... .......... .......... 73% 19.3M 1s
    ##  98100K .......... .......... .......... .......... .......... 73%  102M 1s
    ##  98150K .......... .......... .......... .......... .......... 73%  105M 1s
    ##  98200K .......... .......... .......... .......... .......... 73%  130M 1s
    ##  98250K .......... .......... .......... .......... .......... 73%  123M 1s
    ##  98300K .......... .......... .......... .......... .......... 73% 73.4M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 18.8M 1s
    ##  98400K .......... .......... .......... .......... .......... 73% 86.5M 1s
    ##  98450K .......... .......... .......... .......... .......... 73% 79.8M 1s
    ##  98500K .......... .......... .......... .......... .......... 73%  106M 1s
    ##  98550K .......... .......... .......... .......... .......... 73% 98.1M 1s
    ##  98600K .......... .......... .......... .......... .......... 73% 24.5M 1s
    ##  98650K .......... .......... .......... .......... .......... 73%  123M 1s
    ##  98700K .......... .......... .......... .......... .......... 73% 23.0M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 83.1M 1s
    ##  98800K .......... .......... .......... .......... .......... 73%  128M 1s
    ##  98850K .......... .......... .......... .......... .......... 73% 14.0M 1s
    ##  98900K .......... .......... .......... .......... .......... 73%  115M 1s
    ##  98950K .......... .......... .......... .......... .......... 73% 93.2M 1s
    ##  99000K .......... .......... .......... .......... .......... 73%  131M 1s
    ##  99050K .......... .......... .......... .......... .......... 73%  126M 1s
    ##  99100K .......... .......... .......... .......... .......... 73%  127M 1s
    ##  99150K .......... .......... .......... .......... .......... 73% 23.1M 1s
    ##  99200K .......... .......... .......... .......... .......... 74%  101M 1s
    ##  99250K .......... .......... .......... .......... .......... 74% 51.0M 1s
    ##  99300K .......... .......... .......... .......... .......... 74% 95.2M 1s
    ##  99350K .......... .......... .......... .......... .......... 74% 97.5M 1s
    ##  99400K .......... .......... .......... .......... .......... 74% 35.6M 1s
    ##  99450K .......... .......... .......... .......... .......... 74%  118M 1s
    ##  99500K .......... .......... .......... .......... .......... 74% 27.5M 1s
    ##  99550K .......... .......... .......... .......... .......... 74% 78.7M 1s
    ##  99600K .......... .......... .......... .......... .......... 74% 99.7M 1s
    ##  99650K .......... .......... .......... .......... .......... 74% 55.7M 1s
    ##  99700K .......... .......... .......... .......... .......... 74% 28.2M 1s
    ##  99750K .......... .......... .......... .......... .......... 74% 80.0M 1s
    ##  99800K .......... .......... .......... .......... .......... 74%  107M 1s
    ##  99850K .......... .......... .......... .......... .......... 74%  121M 1s
    ##  99900K .......... .......... .......... .......... .......... 74% 53.9M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 16.3M 1s
    ## 100000K .......... .......... .......... .......... .......... 74%  114M 1s
    ## 100050K .......... .......... .......... .......... .......... 74%  118M 1s
    ## 100100K .......... .......... .......... .......... .......... 74%  137M 1s
    ## 100150K .......... .......... .......... .......... .......... 74%  119M 1s
    ## 100200K .......... .......... .......... .......... .......... 74% 5.86M 1s
    ## 100250K .......... .......... .......... .......... .......... 74%  127M 1s
    ## 100300K .......... .......... .......... .......... .......... 74%  138M 1s
    ## 100350K .......... .......... .......... .......... .......... 74%  115M 1s
    ## 100400K .......... .......... .......... .......... .......... 74%  135M 1s
    ## 100450K .......... .......... .......... .......... .......... 74% 14.4M 1s
    ## 100500K .......... .......... .......... .......... .......... 75% 97.0M 1s
    ## 100550K .......... .......... .......... .......... .......... 75% 93.9M 1s
    ## 100600K .......... .......... .......... .......... .......... 75% 96.5M 1s
    ## 100650K .......... .......... .......... .......... .......... 75%  116M 1s
    ## 100700K .......... .......... .......... .......... .......... 75%  123M 1s
    ## 100750K .......... .......... .......... .......... .......... 75% 21.5M 1s
    ## 100800K .......... .......... .......... .......... .......... 75% 88.4M 1s
    ## 100850K .......... .......... .......... .......... .......... 75% 65.0M 1s
    ## 100900K .......... .......... .......... .......... .......... 75% 95.9M 1s
    ## 100950K .......... .......... .......... .......... .......... 75% 69.2M 1s
    ## 101000K .......... .......... .......... .......... .......... 75% 50.0M 1s
    ## 101050K .......... .......... .......... .......... .......... 75%  102M 1s
    ## 101100K .......... .......... .......... .......... .......... 75% 80.4M 1s
    ## 101150K .......... .......... .......... .......... .......... 75% 21.7M 1s
    ## 101200K .......... .......... .......... .......... .......... 75%  136M 1s
    ## 101250K .......... .......... .......... .......... .......... 75% 67.8M 1s
    ## 101300K .......... .......... .......... .......... .......... 75% 56.7M 1s
    ## 101350K .......... .......... .......... .......... .......... 75% 79.3M 1s
    ## 101400K .......... .......... .......... .......... .......... 75% 33.6M 1s
    ## 101450K .......... .......... .......... .......... .......... 75% 48.8M 1s
    ## 101500K .......... .......... .......... .......... .......... 75%  133M 1s
    ## 101550K .......... .......... .......... .......... .......... 75% 35.9M 1s
    ## 101600K .......... .......... .......... .......... .......... 75%  129M 1s
    ## 101650K .......... .......... .......... .......... .......... 75% 26.0M 1s
    ## 101700K .......... .......... .......... .......... .......... 75%  132M 1s
    ## 101750K .......... .......... .......... .......... .......... 75% 15.6M 1s
    ## 101800K .......... .......... .......... .......... .......... 75% 97.5M 1s
    ## 101850K .......... .......... .......... .......... .......... 76%  116M 1s
    ## 101900K .......... .......... .......... .......... .......... 76% 66.2M 1s
    ## 101950K .......... .......... .......... .......... .......... 76%  100M 1s
    ## 102000K .......... .......... .......... .......... .......... 76% 35.9M 1s
    ## 102050K .......... .......... .......... .......... .......... 76% 19.5M 1s
    ## 102100K .......... .......... .......... .......... .......... 76% 94.5M 1s
    ## 102150K .......... .......... .......... .......... .......... 76%  121M 1s
    ## 102200K .......... .......... .......... .......... .......... 76%  138M 1s
    ## 102250K .......... .......... .......... .......... .......... 76%  134M 1s
    ## 102300K .......... .......... .......... .......... .......... 76%  141M 1s
    ## 102350K .......... .......... .......... .......... .......... 76% 28.6M 1s
    ## 102400K .......... .......... .......... .......... .......... 76% 71.6M 1s
    ## 102450K .......... .......... .......... .......... .......... 76%  109M 1s
    ## 102500K .......... .......... .......... .......... .......... 76%  109M 1s
    ## 102550K .......... .......... .......... .......... .......... 76%  116M 1s
    ## 102600K .......... .......... .......... .......... .......... 76% 9.24M 1s
    ## 102650K .......... .......... .......... .......... .......... 76%  129M 1s
    ## 102700K .......... .......... .......... .......... .......... 76%  120M 1s
    ## 102750K .......... .......... .......... .......... .......... 76%  124M 1s
    ## 102800K .......... .......... .......... .......... .......... 76%  144M 1s
    ## 102850K .......... .......... .......... .......... .......... 76% 23.6M 1s
    ## 102900K .......... .......... .......... .......... .......... 76% 88.8M 1s
    ## 102950K .......... .......... .......... .......... .......... 76% 33.2M 1s
    ## 103000K .......... .......... .......... .......... .......... 76% 96.0M 1s
    ## 103050K .......... .......... .......... .......... .......... 76% 95.0M 1s
    ## 103100K .......... .......... .......... .......... .......... 76% 22.4M 1s
    ## 103150K .......... .......... .......... .......... .......... 76%  118M 1s
    ## 103200K .......... .......... .......... .......... .......... 77% 71.2M 1s
    ## 103250K .......... .......... .......... .......... .......... 77%  139M 1s
    ## 103300K .......... .......... .......... .......... .......... 77% 78.4M 1s
    ## 103350K .......... .......... .......... .......... .......... 77% 15.4M 1s
    ## 103400K .......... .......... .......... .......... .......... 77%  115M 1s
    ## 103450K .......... .......... .......... .......... .......... 77%  142M 1s
    ## 103500K .......... .......... .......... .......... .......... 77%  113M 1s
    ## 103550K .......... .......... .......... .......... .......... 77%  116M 1s
    ## 103600K .......... .......... .......... .......... .......... 77% 4.45M 1s
    ## 103650K .......... .......... .......... .......... .......... 77% 50.7M 1s
    ## 103700K .......... .......... .......... .......... .......... 77% 49.9M 1s
    ## 103750K .......... .......... .......... .......... .......... 77% 93.7M 1s
    ## 103800K .......... .......... .......... .......... .......... 77%  101M 1s
    ## 103850K .......... .......... .......... .......... .......... 77%  117M 1s
    ## 103900K .......... .......... .......... .......... .......... 77% 11.6M 1s
    ## 103950K .......... .......... .......... .......... .......... 77% 87.0M 1s
    ## 104000K .......... .......... .......... .......... .......... 77%  142M 1s
    ## 104050K .......... .......... .......... .......... .......... 77% 97.1M 1s
    ## 104100K .......... .......... .......... .......... .......... 77%  124M 1s
    ## 104150K .......... .......... .......... .......... .......... 77%  123M 1s
    ## 104200K .......... .......... .......... .......... .......... 77% 21.8M 1s
    ## 104250K .......... .......... .......... .......... .......... 77% 79.1M 1s
    ## 104300K .......... .......... .......... .......... .......... 77% 71.2M 1s
    ## 104350K .......... .......... .......... .......... .......... 77% 55.4M 1s
    ## 104400K .......... .......... .......... .......... .......... 77% 33.2M 1s
    ## 104450K .......... .......... .......... .......... .......... 77%  115M 1s
    ## 104500K .......... .......... .......... .......... .......... 77% 22.6M 1s
    ## 104550K .......... .......... .......... .......... .......... 78% 91.0M 1s
    ## 104600K .......... .......... .......... .......... .......... 78%  132M 1s
    ## 104650K .......... .......... .......... .......... .......... 78%  109M 1s
    ## 104700K .......... .......... .......... .......... .......... 78%  131M 1s
    ## 104750K .......... .......... .......... .......... .......... 78% 41.7M 1s
    ## 104800K .......... .......... .......... .......... .......... 78%  106M 1s
    ## 104850K .......... .......... .......... .......... .......... 78%  117M 1s
    ## 104900K .......... .......... .......... .......... .......... 78% 91.2M 1s
    ## 104950K .......... .......... .......... .......... .......... 78%  105M 1s
    ## 105000K .......... .......... .......... .......... .......... 78% 35.8M 1s
    ## 105050K .......... .......... .......... .......... .......... 78%  113M 1s
    ## 105100K .......... .......... .......... .......... .......... 78%  121M 1s
    ## 105150K .......... .......... .......... .......... .......... 78% 10.4M 1s
    ## 105200K .......... .......... .......... .......... .......... 78%  130M 1s
    ## 105250K .......... .......... .......... .......... .......... 78%  139M 1s
    ## 105300K .......... .......... .......... .......... .......... 78%  129M 1s
    ## 105350K .......... .......... .......... .......... .......... 78%  128M 1s
    ## 105400K .......... .......... .......... .......... .......... 78% 24.5M 1s
    ## 105450K .......... .......... .......... .......... .......... 78% 35.1M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 48.7M 1s
    ## 105550K .......... .......... .......... .......... .......... 78% 89.0M 1s
    ## 105600K .......... .......... .......... .......... .......... 78% 61.9M 1s
    ## 105650K .......... .......... .......... .......... .......... 78% 49.5M 1s
    ## 105700K .......... .......... .......... .......... .......... 78%  122M 1s
    ## 105750K .......... .......... .......... .......... .......... 78% 24.5M 1s
    ## 105800K .......... .......... .......... .......... .......... 78% 99.6M 1s
    ## 105850K .......... .......... .......... .......... .......... 78%  127M 1s
    ## 105900K .......... .......... .......... .......... .......... 79% 31.8M 1s
    ## 105950K .......... .......... .......... .......... .......... 79%  104M 1s
    ## 106000K .......... .......... .......... .......... .......... 79% 38.0M 1s
    ## 106050K .......... .......... .......... .......... .......... 79%  111M 1s
    ## 106100K .......... .......... .......... .......... .......... 79% 30.8M 1s
    ## 106150K .......... .......... .......... .......... .......... 79% 95.0M 1s
    ## 106200K .......... .......... .......... .......... .......... 79%  122M 1s
    ## 106250K .......... .......... .......... .......... .......... 79%  113M 1s
    ## 106300K .......... .......... .......... .......... .......... 79%  105M 1s
    ## 106350K .......... .......... .......... .......... .......... 79% 18.6M 1s
    ## 106400K .......... .......... .......... .......... .......... 79% 27.4M 1s
    ## 106450K .......... .......... .......... .......... .......... 79%  125M 1s
    ## 106500K .......... .......... .......... .......... .......... 79%  151M 1s
    ## 106550K .......... .......... .......... .......... .......... 79%  120M 1s
    ## 106600K .......... .......... .......... .......... .......... 79%  138M 1s
    ## 106650K .......... .......... .......... .......... .......... 79% 39.9M 1s
    ## 106700K .......... .......... .......... .......... .......... 79%  110M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 37.3M 1s
    ## 106800K .......... .......... .......... .......... .......... 79%  115M 1s
    ## 106850K .......... .......... .......... .......... .......... 79% 53.5M 1s
    ## 106900K .......... .......... .......... .......... .......... 79% 40.3M 1s
    ## 106950K .......... .......... .......... .......... .......... 79% 78.3M 1s
    ## 107000K .......... .......... .......... .......... .......... 79% 90.1M 1s
    ## 107050K .......... .......... .......... .......... .......... 79% 39.3M 1s
    ## 107100K .......... .......... .......... .......... .......... 79% 71.4M 1s
    ## 107150K .......... .......... .......... .......... .......... 79% 46.5M 1s
    ## 107200K .......... .......... .......... .......... .......... 79% 60.6M 1s
    ## 107250K .......... .......... .......... .......... .......... 80%  104M 1s
    ## 107300K .......... .......... .......... .......... .......... 80% 47.9M 1s
    ## 107350K .......... .......... .......... .......... .......... 80% 49.7M 1s
    ## 107400K .......... .......... .......... .......... .......... 80% 55.7M 1s
    ## 107450K .......... .......... .......... .......... .......... 80%  105M 1s
    ## 107500K .......... .......... .......... .......... .......... 80% 65.6M 1s
    ## 107550K .......... .......... .......... .......... .......... 80% 55.4M 1s
    ## 107600K .......... .......... .......... .......... .......... 80% 58.4M 1s
    ## 107650K .......... .......... .......... .......... .......... 80% 51.9M 1s
    ## 107700K .......... .......... .......... .......... .......... 80%  112M 1s
    ## 107750K .......... .......... .......... .......... .......... 80% 66.2M 1s
    ## 107800K .......... .......... .......... .......... .......... 80% 45.3M 1s
    ## 107850K .......... .......... .......... .......... .......... 80% 76.6M 1s
    ## 107900K .......... .......... .......... .......... .......... 80%  119M 1s
    ## 107950K .......... .......... .......... .......... .......... 80% 45.5M 1s
    ## 108000K .......... .......... .......... .......... .......... 80% 56.5M 1s
    ## 108050K .......... .......... .......... .......... .......... 80% 55.6M 1s
    ## 108100K .......... .......... .......... .......... .......... 80% 57.8M 1s
    ## 108150K .......... .......... .......... .......... .......... 80% 61.1M 1s
    ## 108200K .......... .......... .......... .......... .......... 80%  106M 1s
    ## 108250K .......... .......... .......... .......... .......... 80% 78.4M 1s
    ## 108300K .......... .......... .......... .......... .......... 80% 51.5M 1s
    ## 108350K .......... .......... .......... .......... .......... 80% 65.8M 1s
    ## 108400K .......... .......... .......... .......... .......... 80% 56.5M 1s
    ## 108450K .......... .......... .......... .......... .......... 80%  111M 1s
    ## 108500K .......... .......... .......... .......... .......... 80% 62.0M 1s
    ## 108550K .......... .......... .......... .......... .......... 81% 64.8M 1s
    ## 108600K .......... .......... .......... .......... .......... 81% 55.7M 1s
    ## 108650K .......... .......... .......... .......... .......... 81%  115M 1s
    ## 108700K .......... .......... .......... .......... .......... 81% 56.5M 1s
    ## 108750K .......... .......... .......... .......... .......... 81% 62.4M 1s
    ## 108800K .......... .......... .......... .......... .......... 81% 55.2M 1s
    ## 108850K .......... .......... .......... .......... .......... 81% 77.0M 1s
    ## 108900K .......... .......... .......... .......... .......... 81%  111M 1s
    ## 108950K .......... .......... .......... .......... .......... 81% 57.3M 1s
    ## 109000K .......... .......... .......... .......... .......... 81% 66.8M 1s
    ## 109050K .......... .......... .......... .......... .......... 81% 84.6M 1s
    ## 109100K .......... .......... .......... .......... .......... 81%  111M 1s
    ## 109150K .......... .......... .......... .......... .......... 81% 56.3M 1s
    ## 109200K .......... .......... .......... .......... .......... 81% 66.7M 1s
    ## 109250K .......... .......... .......... .......... .......... 81% 51.2M 1s
    ## 109300K .......... .......... .......... .......... .......... 81% 93.2M 1s
    ## 109350K .......... .......... .......... .......... .......... 81%  111M 1s
    ## 109400K .......... .......... .......... .......... .......... 81% 43.4M 1s
    ## 109450K .......... .......... .......... .......... .......... 81% 68.0M 1s
    ## 109500K .......... .......... .......... .......... .......... 81% 50.8M 1s
    ## 109550K .......... .......... .......... .......... .......... 81% 58.1M 1s
    ## 109600K .......... .......... .......... .......... .......... 81%  111M 1s
    ## 109650K .......... .......... .......... .......... .......... 81% 65.1M 1s
    ## 109700K .......... .......... .......... .......... .......... 81% 73.9M 1s
    ## 109750K .......... .......... .......... .......... .......... 81% 69.4M 1s
    ## 109800K .......... .......... .......... .......... .......... 81% 77.8M 1s
    ## 109850K .......... .......... .......... .......... .......... 81% 43.9M 1s
    ## 109900K .......... .......... .......... .......... .......... 82% 53.8M 1s
    ## 109950K .......... .......... .......... .......... .......... 82% 48.5M 1s
    ## 110000K .......... .......... .......... .......... .......... 82% 59.9M 1s
    ## 110050K .......... .......... .......... .......... .......... 82% 77.7M 1s
    ## 110100K .......... .......... .......... .......... .......... 82% 56.3M 1s
    ## 110150K .......... .......... .......... .......... .......... 82% 47.6M 1s
    ## 110200K .......... .......... .......... .......... .......... 82% 56.7M 1s
    ## 110250K .......... .......... .......... .......... .......... 82% 57.5M 1s
    ## 110300K .......... .......... .......... .......... .......... 82% 84.0M 1s
    ## 110350K .......... .......... .......... .......... .......... 82% 52.0M 1s
    ## 110400K .......... .......... .......... .......... .......... 82% 52.6M 1s
    ## 110450K .......... .......... .......... .......... .......... 82% 55.2M 1s
    ## 110500K .......... .......... .......... .......... .......... 82% 76.1M 1s
    ## 110550K .......... .......... .......... .......... .......... 82% 87.3M 1s
    ## 110600K .......... .......... .......... .......... .......... 82% 50.4M 1s
    ## 110650K .......... .......... .......... .......... .......... 82% 71.3M 1s
    ## 110700K .......... .......... .......... .......... .......... 82% 59.1M 1s
    ## 110750K .......... .......... .......... .......... .......... 82% 51.1M 1s
    ## 110800K .......... .......... .......... .......... .......... 82% 95.9M 1s
    ## 110850K .......... .......... .......... .......... .......... 82% 52.7M 1s
    ## 110900K .......... .......... .......... .......... .......... 82% 60.8M 1s
    ## 110950K .......... .......... .......... .......... .......... 82% 66.9M 1s
    ## 111000K .......... .......... .......... .......... .......... 82% 86.6M 1s
    ## 111050K .......... .......... .......... .......... .......... 82% 66.5M 1s
    ## 111100K .......... .......... .......... .......... .......... 82% 52.9M 1s
    ## 111150K .......... .......... .......... .......... .......... 82% 60.8M 1s
    ## 111200K .......... .......... .......... .......... .......... 82% 61.1M 1s
    ## 111250K .......... .......... .......... .......... .......... 83%  102M 1s
    ## 111300K .......... .......... .......... .......... .......... 83% 46.3M 1s
    ## 111350K .......... .......... .......... .......... .......... 83% 57.8M 1s
    ## 111400K .......... .......... .......... .......... .......... 83% 57.0M 1s
    ## 111450K .......... .......... .......... .......... .......... 83% 90.7M 1s
    ## 111500K .......... .......... .......... .......... .......... 83% 68.8M 1s
    ## 111550K .......... .......... .......... .......... .......... 83% 56.3M 1s
    ## 111600K .......... .......... .......... .......... .......... 83% 25.5M 1s
    ## 111650K .......... .......... .......... .......... .......... 83%  104M 1s
    ## 111700K .......... .......... .......... .......... .......... 83%  105M 1s
    ## 111750K .......... .......... .......... .......... .......... 83%  122M 1s
    ## 111800K .......... .......... .......... .......... .......... 83%  121M 1s
    ## 111850K .......... .......... .......... .......... .......... 83% 27.8M 1s
    ## 111900K .......... .......... .......... .......... .......... 83%  112M 1s
    ## 111950K .......... .......... .......... .......... .......... 83% 58.0M 1s
    ## 112000K .......... .......... .......... .......... .......... 83%  122M 1s
    ## 112050K .......... .......... .......... .......... .......... 83% 76.5M 1s
    ## 112100K .......... .......... .......... .......... .......... 83% 27.6M 1s
    ## 112150K .......... .......... .......... .......... .......... 83% 84.9M 1s
    ## 112200K .......... .......... .......... .......... .......... 83% 98.3M 1s
    ## 112250K .......... .......... .......... .......... .......... 83% 52.2M 1s
    ## 112300K .......... .......... .......... .......... .......... 83% 62.2M 1s
    ## 112350K .......... .......... .......... .......... .......... 83% 39.9M 1s
    ## 112400K .......... .......... .......... .......... .......... 83% 70.5M 1s
    ## 112450K .......... .......... .......... .......... .......... 83% 81.7M 1s
    ## 112500K .......... .......... .......... .......... .......... 83% 22.6M 1s
    ## 112550K .......... .......... .......... .......... .......... 83%  117M 1s
    ## 112600K .......... .......... .......... .......... .......... 84% 49.4M 1s
    ## 112650K .......... .......... .......... .......... .......... 84% 89.8M 1s
    ## 112700K .......... .......... .......... .......... .......... 84% 40.7M 1s
    ## 112750K .......... .......... .......... .......... .......... 84% 76.7M 1s
    ## 112800K .......... .......... .......... .......... .......... 84% 18.8M 1s
    ## 112850K .......... .......... .......... .......... .......... 84%  139M 1s
    ## 112900K .......... .......... .......... .......... .......... 84%  150M 1s
    ## 112950K .......... .......... .......... .......... .......... 84% 28.2M 1s
    ## 113000K .......... .......... .......... .......... .......... 84%  149M 1s
    ## 113050K .......... .......... .......... .......... .......... 84%  132M 1s
    ## 113100K .......... .......... .......... .......... .......... 84% 28.3M 1s
    ## 113150K .......... .......... .......... .......... .......... 84% 89.6M 1s
    ## 113200K .......... .......... .......... .......... .......... 84% 71.5M 1s
    ## 113250K .......... .......... .......... .......... .......... 84% 91.7M 1s
    ## 113300K .......... .......... .......... .......... .......... 84% 84.6M 1s
    ## 113350K .......... .......... .......... .......... .......... 84% 95.5M 1s
    ## 113400K .......... .......... .......... .......... .......... 84% 22.4M 1s
    ## 113450K .......... .......... .......... .......... .......... 84%  127M 1s
    ## 113500K .......... .......... .......... .......... .......... 84% 36.7M 1s
    ## 113550K .......... .......... .......... .......... .......... 84% 78.1M 1s
    ## 113600K .......... .......... .......... .......... .......... 84%  141M 1s
    ## 113650K .......... .......... .......... .......... .......... 84% 40.4M 1s
    ## 113700K .......... .......... .......... .......... .......... 84%  105M 1s
    ## 113750K .......... .......... .......... .......... .......... 84% 43.6M 1s
    ## 113800K .......... .......... .......... .......... .......... 84% 90.4M 1s
    ## 113850K .......... .......... .......... .......... .......... 84%  132M 1s
    ## 113900K .......... .......... .......... .......... .......... 84% 26.8M 1s
    ## 113950K .......... .......... .......... .......... .......... 85% 83.9M 1s
    ## 114000K .......... .......... .......... .......... .......... 85% 35.0M 1s
    ## 114050K .......... .......... .......... .......... .......... 85% 97.2M 1s
    ## 114100K .......... .......... .......... .......... .......... 85% 64.0M 1s
    ## 114150K .......... .......... .......... .......... .......... 85% 63.1M 1s
    ## 114200K .......... .......... .......... .......... .......... 85% 57.9M 1s
    ## 114250K .......... .......... .......... .......... .......... 85%  120M 1s
    ## 114300K .......... .......... .......... .......... .......... 85% 45.2M 1s
    ## 114350K .......... .......... .......... .......... .......... 85% 42.8M 1s
    ## 114400K .......... .......... .......... .......... .......... 85%  131M 1s
    ## 114450K .......... .......... .......... .......... .......... 85%  102M 1s
    ## 114500K .......... .......... .......... .......... .......... 85%  140M 1s
    ## 114550K .......... .......... .......... .......... .......... 85% 34.6M 1s
    ## 114600K .......... .......... .......... .......... .......... 85% 27.4M 1s
    ## 114650K .......... .......... .......... .......... .......... 85%  108M 1s
    ## 114700K .......... .......... .......... .......... .......... 85%  109M 1s
    ## 114750K .......... .......... .......... .......... .......... 85%  101M 1s
    ## 114800K .......... .......... .......... .......... .......... 85%  111M 1s
    ## 114850K .......... .......... .......... .......... .......... 85% 25.6M 1s
    ## 114900K .......... .......... .......... .......... .......... 85% 55.7M 1s
    ## 114950K .......... .......... .......... .......... .......... 85% 68.7M 1s
    ## 115000K .......... .......... .......... .......... .......... 85%  107M 1s
    ## 115050K .......... .......... .......... .......... .......... 85% 85.9M 1s
    ## 115100K .......... .......... .......... .......... .......... 85% 57.5M 1s
    ## 115150K .......... .......... .......... .......... .......... 85% 46.8M 1s
    ## 115200K .......... .......... .......... .......... .......... 85% 35.3M 1s
    ## 115250K .......... .......... .......... .......... .......... 86% 87.0M 1s
    ## 115300K .......... .......... .......... .......... .......... 86%  137M 1s
    ## 115350K .......... .......... .......... .......... .......... 86% 40.2M 1s
    ## 115400K .......... .......... .......... .......... .......... 86%  128M 1s
    ## 115450K .......... .......... .......... .......... .......... 86%  128M 1s
    ## 115500K .......... .......... .......... .......... .......... 86% 13.7M 1s
    ## 115550K .......... .......... .......... .......... .......... 86%  120M 1s
    ## 115600K .......... .......... .......... .......... .......... 86%  148M 1s
    ## 115650K .......... .......... .......... .......... .......... 86%  127M 1s
    ## 115700K .......... .......... .......... .......... .......... 86%  146M 1s
    ## 115750K .......... .......... .......... .......... .......... 86% 17.7M 1s
    ## 115800K .......... .......... .......... .......... .......... 86%  128M 1s
    ## 115850K .......... .......... .......... .......... .......... 86% 62.8M 1s
    ## 115900K .......... .......... .......... .......... .......... 86% 95.8M 1s
    ## 115950K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 116000K .......... .......... .......... .......... .......... 86% 28.6M 1s
    ## 116050K .......... .......... .......... .......... .......... 86% 94.6M 1s
    ## 116100K .......... .......... .......... .......... .......... 86% 39.1M 1s
    ## 116150K .......... .......... .......... .......... .......... 86% 70.5M 1s
    ## 116200K .......... .......... .......... .......... .......... 86%  120M 1s
    ## 116250K .......... .......... .......... .......... .......... 86% 8.86M 1s
    ## 116300K .......... .......... .......... .......... .......... 86%  105M 1s
    ## 116350K .......... .......... .......... .......... .......... 86% 79.7M 1s
    ## 116400K .......... .......... .......... .......... .......... 86%  110M 1s
    ## 116450K .......... .......... .......... .......... .......... 86%  136M 1s
    ## 116500K .......... .......... .......... .......... .......... 86%  140M 1s
    ## 116550K .......... .......... .......... .......... .......... 86% 35.8M 1s
    ## 116600K .......... .......... .......... .......... .......... 87% 35.4M 1s
    ## 116650K .......... .......... .......... .......... .......... 87% 49.0M 1s
    ## 116700K .......... .......... .......... .......... .......... 87% 78.7M 0s
    ## 116750K .......... .......... .......... .......... .......... 87% 85.2M 0s
    ## 116800K .......... .......... .......... .......... .......... 87%  118M 0s
    ## 116850K .......... .......... .......... .......... .......... 87% 11.9M 0s
    ## 116900K .......... .......... .......... .......... .......... 87% 89.6M 0s
    ## 116950K .......... .......... .......... .......... .......... 87%  108M 0s
    ## 117000K .......... .......... .......... .......... .......... 87%  139M 0s
    ## 117050K .......... .......... .......... .......... .......... 87%  137M 0s
    ## 117100K .......... .......... .......... .......... .......... 87%  136M 0s
    ## 117150K .......... .......... .......... .......... .......... 87% 8.56M 0s
    ## 117200K .......... .......... .......... .......... .......... 87%  123M 0s
    ## 117250K .......... .......... .......... .......... .......... 87% 45.5M 0s
    ## 117300K .......... .......... .......... .......... .......... 87%  127M 0s
    ## 117350K .......... .......... .......... .......... .......... 87%  104M 0s
    ## 117400K .......... .......... .......... .......... .......... 87%  136M 0s
    ## 117450K .......... .......... .......... .......... .......... 87%  117M 0s
    ## 117500K .......... .......... .......... .......... .......... 87%  103M 0s
    ## 117550K .......... .......... .......... .......... .......... 87% 63.0M 0s
    ## 117600K .......... .......... .......... .......... .......... 87% 49.0M 0s
    ## 117650K .......... .......... .......... .......... .......... 87%  104M 0s
    ## 117700K .......... .......... .......... .......... .......... 87% 23.8M 0s
    ## 117750K .......... .......... .......... .......... .......... 87%  104M 0s
    ## 117800K .......... .......... .......... .......... .......... 87%  135M 0s
    ## 117850K .......... .......... .......... .......... .......... 87% 17.0M 0s
    ## 117900K .......... .......... .......... .......... .......... 87%  137M 0s
    ## 117950K .......... .......... .......... .......... .......... 88%  108M 0s
    ## 118000K .......... .......... .......... .......... .......... 88% 35.2M 0s
    ## 118050K .......... .......... .......... .......... .......... 88% 86.2M 0s
    ## 118100K .......... .......... .......... .......... .......... 88% 36.1M 0s
    ## 118150K .......... .......... .......... .......... .......... 88% 93.2M 0s
    ## 118200K .......... .......... .......... .......... .......... 88%  116M 0s
    ## 118250K .......... .......... .......... .......... .......... 88%  138M 0s
    ## 118300K .......... .......... .......... .......... .......... 88% 67.3M 0s
    ## 118350K .......... .......... .......... .......... .......... 88% 99.1M 0s
    ## 118400K .......... .......... .......... .......... .......... 88% 32.2M 0s
    ## 118450K .......... .......... .......... .......... .......... 88% 89.1M 0s
    ## 118500K .......... .......... .......... .......... .......... 88% 95.0M 0s
    ## 118550K .......... .......... .......... .......... .......... 88% 73.4M 0s
    ## 118600K .......... .......... .......... .......... .......... 88% 95.4M 0s
    ## 118650K .......... .......... .......... .......... .......... 88% 31.3M 0s
    ## 118700K .......... .......... .......... .......... .......... 88% 95.9M 0s
    ## 118750K .......... .......... .......... .......... .......... 88% 39.1M 0s
    ## 118800K .......... .......... .......... .......... .......... 88% 93.4M 0s
    ## 118850K .......... .......... .......... .......... .......... 88% 58.9M 0s
    ## 118900K .......... .......... .......... .......... .......... 88% 87.5M 0s
    ## 118950K .......... .......... .......... .......... .......... 88% 77.1M 0s
    ## 119000K .......... .......... .......... .......... .......... 88% 98.5M 0s
    ## 119050K .......... .......... .......... .......... .......... 88% 26.1M 0s
    ## 119100K .......... .......... .......... .......... .......... 88%  137M 0s
    ## 119150K .......... .......... .......... .......... .......... 88% 16.6M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 35.9M 0s
    ## 119250K .......... .......... .......... .......... .......... 88% 54.1M 0s
    ## 119300K .......... .......... .......... .......... .......... 89% 97.9M 0s
    ## 119350K .......... .......... .......... .......... .......... 89%  105M 0s
    ## 119400K .......... .......... .......... .......... .......... 89% 38.6M 0s
    ## 119450K .......... .......... .......... .......... .......... 89%  126M 0s
    ## 119500K .......... .......... .......... .......... .......... 89%  142M 0s
    ## 119550K .......... .......... .......... .......... .......... 89% 47.7M 0s
    ## 119600K .......... .......... .......... .......... .......... 89% 99.2M 0s
    ## 119650K .......... .......... .......... .......... .......... 89%  108M 0s
    ## 119700K .......... .......... .......... .......... .......... 89%  137M 0s
    ## 119750K .......... .......... .......... .......... .......... 89% 45.9M 0s
    ## 119800K .......... .......... .......... .......... .......... 89% 28.1M 0s
    ## 119850K .......... .......... .......... .......... .......... 89% 99.2M 0s
    ## 119900K .......... .......... .......... .......... .......... 89%  128M 0s
    ## 119950K .......... .......... .......... .......... .......... 89%  109M 0s
    ## 120000K .......... .......... .......... .......... .......... 89%  112M 0s
    ## 120050K .......... .......... .......... .......... .......... 89% 28.1M 0s
    ## 120100K .......... .......... .......... .......... .......... 89%  143M 0s
    ## 120150K .......... .......... .......... .......... .......... 89% 96.5M 0s
    ## 120200K .......... .......... .......... .......... .......... 89%  109M 0s
    ## 120250K .......... .......... .......... .......... .......... 89%  133M 0s
    ## 120300K .......... .......... .......... .......... .......... 89% 13.0M 0s
    ## 120350K .......... .......... .......... .......... .......... 89% 97.3M 0s
    ## 120400K .......... .......... .......... .......... .......... 89%  149M 0s
    ## 120450K .......... .......... .......... .......... .......... 89%  147M 0s
    ## 120500K .......... .......... .......... .......... .......... 89%  142M 0s
    ## 120550K .......... .......... .......... .......... .......... 89% 27.3M 0s
    ## 120600K .......... .......... .......... .......... .......... 89%  113M 0s
    ## 120650K .......... .......... .......... .......... .......... 90% 42.8M 0s
    ## 120700K .......... .......... .......... .......... .......... 90% 99.1M 0s
    ## 120750K .......... .......... .......... .......... .......... 90% 87.3M 0s
    ## 120800K .......... .......... .......... .......... .......... 90%  120M 0s
    ## 120850K .......... .......... .......... .......... .......... 90% 43.3M 0s
    ## 120900K .......... .......... .......... .......... .......... 90%  126M 0s
    ## 120950K .......... .......... .......... .......... .......... 90% 39.0M 0s
    ## 121000K .......... .......... .......... .......... .......... 90%  124M 0s
    ## 121050K .......... .......... .......... .......... .......... 90% 72.0M 0s
    ## 121100K .......... .......... .......... .......... .......... 90%  108M 0s
    ## 121150K .......... .......... .......... .......... .......... 90% 41.1M 0s
    ## 121200K .......... .......... .......... .......... .......... 90% 37.4M 0s
    ## 121250K .......... .......... .......... .......... .......... 90%  102M 0s
    ## 121300K .......... .......... .......... .......... .......... 90% 29.0M 0s
    ## 121350K .......... .......... .......... .......... .......... 90%  122M 0s
    ## 121400K .......... .......... .......... .......... .......... 90% 35.3M 0s
    ## 121450K .......... .......... .......... .......... .......... 90%  139M 0s
    ## 121500K .......... .......... .......... .......... .......... 90% 45.9M 0s
    ## 121550K .......... .......... .......... .......... .......... 90%  123M 0s
    ## 121600K .......... .......... .......... .......... .......... 90%  154M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 19.3M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 95.1M 0s
    ## 121750K .......... .......... .......... .......... .......... 90% 27.4M 0s
    ## 121800K .......... .......... .......... .......... .......... 90% 71.2M 0s
    ## 121850K .......... .......... .......... .......... .......... 90%  109M 0s
    ## 121900K .......... .......... .......... .......... .......... 90%  139M 0s
    ## 121950K .......... .......... .......... .......... .......... 91% 99.6M 0s
    ## 122000K .......... .......... .......... .......... .......... 91% 70.5M 0s
    ## 122050K .......... .......... .......... .......... .......... 91%  105M 0s
    ## 122100K .......... .......... .......... .......... .......... 91% 56.8M 0s
    ## 122150K .......... .......... .......... .......... .......... 91% 38.0M 0s
    ## 122200K .......... .......... .......... .......... .......... 91%  112M 0s
    ## 122250K .......... .......... .......... .......... .......... 91% 54.9M 0s
    ## 122300K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 122350K .......... .......... .......... .......... .......... 91% 63.1M 0s
    ## 122400K .......... .......... .......... .......... .......... 91% 39.2M 0s
    ## 122450K .......... .......... .......... .......... .......... 91%  134M 0s
    ## 122500K .......... .......... .......... .......... .......... 91%  146M 0s
    ## 122550K .......... .......... .......... .......... .......... 91% 38.1M 0s
    ## 122600K .......... .......... .......... .......... .......... 91% 80.7M 0s
    ## 122650K .......... .......... .......... .......... .......... 91%  100M 0s
    ## 122700K .......... .......... .......... .......... .......... 91% 35.6M 0s
    ## 122750K .......... .......... .......... .......... .......... 91% 31.0M 0s
    ## 122800K .......... .......... .......... .......... .......... 91% 31.8M 0s
    ## 122850K .......... .......... .......... .......... .......... 91%  127M 0s
    ## 122900K .......... .......... .......... .......... .......... 91%  138M 0s
    ## 122950K .......... .......... .......... .......... .......... 91%  133M 0s
    ## 123000K .......... .......... .......... .......... .......... 91%  152M 0s
    ## 123050K .......... .......... .......... .......... .......... 91% 23.4M 0s
    ## 123100K .......... .......... .......... .......... .......... 91%  109M 0s
    ## 123150K .......... .......... .......... .......... .......... 91%  124M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 70.6M 0s
    ## 123250K .......... .......... .......... .......... .......... 91% 81.1M 0s
    ## 123300K .......... .......... .......... .......... .......... 92% 26.2M 0s
    ## 123350K .......... .......... .......... .......... .......... 92% 71.4M 0s
    ## 123400K .......... .......... .......... .......... .......... 92% 76.1M 0s
    ## 123450K .......... .......... .......... .......... .......... 92%  102M 0s
    ## 123500K .......... .......... .......... .......... .......... 92% 87.4M 0s
    ## 123550K .......... .......... .......... .......... .......... 92% 95.0M 0s
    ## 123600K .......... .......... .......... .......... .......... 92% 23.7M 0s
    ## 123650K .......... .......... .......... .......... .......... 92%  105M 0s
    ## 123700K .......... .......... .......... .......... .......... 92%  106M 0s
    ## 123750K .......... .......... .......... .......... .......... 92% 69.9M 0s
    ## 123800K .......... .......... .......... .......... .......... 92%  110M 0s
    ## 123850K .......... .......... .......... .......... .......... 92%  119M 0s
    ## 123900K .......... .......... .......... .......... .......... 92%  116M 0s
    ## 123950K .......... .......... .......... .......... .......... 92% 16.5M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 82.6M 0s
    ## 124050K .......... .......... .......... .......... .......... 92%  117M 0s
    ## 124100K .......... .......... .......... .......... .......... 92% 85.8M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 75.9M 0s
    ## 124200K .......... .......... .......... .......... .......... 92%  114M 0s
    ## 124250K .......... .......... .......... .......... .......... 92%  108M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 75.9M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 48.0M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 76.2M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 79.7M 0s
    ## 124500K .......... .......... .......... .......... .......... 92%  110M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 34.2M 0s
    ## 124600K .......... .......... .......... .......... .......... 92% 20.3M 0s
    ## 124650K .......... .......... .......... .......... .......... 93%  108M 0s
    ## 124700K .......... .......... .......... .......... .......... 93%  118M 0s
    ## 124750K .......... .......... .......... .......... .......... 93%  100M 0s
    ## 124800K .......... .......... .......... .......... .......... 93% 97.9M 0s
    ## 124850K .......... .......... .......... .......... .......... 93%  116M 0s
    ## 124900K .......... .......... .......... .......... .......... 93%  121M 0s
    ## 124950K .......... .......... .......... .......... .......... 93%  100M 0s
    ## 125000K .......... .......... .......... .......... .......... 93% 58.6M 0s
    ## 125050K .......... .......... .......... .......... .......... 93% 80.7M 0s
    ## 125100K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125150K .......... .......... .......... .......... .......... 93% 83.7M 0s
    ## 125200K .......... .......... .......... .......... .......... 93% 97.5M 0s
    ## 125250K .......... .......... .......... .......... .......... 93%  121M 0s
    ## 125300K .......... .......... .......... .......... .......... 93% 83.4M 0s
    ## 125350K .......... .......... .......... .......... .......... 93% 86.8M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 88.5M 0s
    ## 125450K .......... .......... .......... .......... .......... 93%  107M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 95.6M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 86.6M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 95.4M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 32.4M 0s
    ## 125700K .......... .......... .......... .......... .......... 93% 84.0M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 61.2M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 54.8M 0s
    ## 125850K .......... .......... .......... .......... .......... 93% 31.9M 0s
    ## 125900K .......... .......... .......... .......... .......... 93%  114M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 92.5M 0s
    ## 126000K .......... .......... .......... .......... .......... 94%  117M 0s
    ## 126050K .......... .......... .......... .......... .......... 94%  128M 0s
    ## 126100K .......... .......... .......... .......... .......... 94%  127M 0s
    ## 126150K .......... .......... .......... .......... .......... 94% 95.6M 0s
    ## 126200K .......... .......... .......... .......... .......... 94% 35.2M 0s
    ## 126250K .......... .......... .......... .......... .......... 94% 99.4M 0s
    ## 126300K .......... .......... .......... .......... .......... 94%  108M 0s
    ## 126350K .......... .......... .......... .......... .......... 94% 83.9M 0s
    ## 126400K .......... .......... .......... .......... .......... 94%  116M 0s
    ## 126450K .......... .......... .......... .......... .......... 94% 28.0M 0s
    ## 126500K .......... .......... .......... .......... .......... 94% 80.9M 0s
    ## 126550K .......... .......... .......... .......... .......... 94% 95.7M 0s
    ## 126600K .......... .......... .......... .......... .......... 94% 93.6M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 96.9M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  128M 0s
    ## 126750K .......... .......... .......... .......... .......... 94%  107M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 69.7M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 74.4M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 96.3M 0s
    ## 126950K .......... .......... .......... .......... .......... 94%  103M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 37.3M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 86.9M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 37.8M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 64.9M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 75.2M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 94.3M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 88.6M 0s
    ## 127350K .......... .......... .......... .......... .......... 95% 82.1M 0s
    ## 127400K .......... .......... .......... .......... .......... 95% 94.6M 0s
    ## 127450K .......... .......... .......... .......... .......... 95% 97.7M 0s
    ## 127500K .......... .......... .......... .......... .......... 95% 73.6M 0s
    ## 127550K .......... .......... .......... .......... .......... 95% 66.4M 0s
    ## 127600K .......... .......... .......... .......... .......... 95% 93.9M 0s
    ## 127650K .......... .......... .......... .......... .......... 95% 76.9M 0s
    ## 127700K .......... .......... .......... .......... .......... 95% 91.6M 0s
    ## 127750K .......... .......... .......... .......... .......... 95% 76.9M 0s
    ## 127800K .......... .......... .......... .......... .......... 95% 76.8M 0s
    ## 127850K .......... .......... .......... .......... .......... 95%  101M 0s
    ## 127900K .......... .......... .......... .......... .......... 95% 95.6M 0s
    ## 127950K .......... .......... .......... .......... .......... 95% 58.4M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 92.9M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 76.3M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 89.7M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 87.1M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 94.1M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 77.9M 0s
    ## 128300K .......... .......... .......... .......... .......... 95%  102M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 65.0M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 82.8M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 91.8M 0s
    ## 128500K .......... .......... .......... .......... .......... 95%  111M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 74.9M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 74.6M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 98.4M 0s
    ## 128700K .......... .......... .......... .......... .......... 96%  103M 0s
    ## 128750K .......... .......... .......... .......... .......... 96% 61.0M 0s
    ## 128800K .......... .......... .......... .......... .......... 96%  116M 0s
    ## 128850K .......... .......... .......... .......... .......... 96%  102M 0s
    ## 128900K .......... .......... .......... .......... .......... 96% 56.3M 0s
    ## 128950K .......... .......... .......... .......... .......... 96% 30.6M 0s
    ## 129000K .......... .......... .......... .......... .......... 96%  119M 0s
    ## 129050K .......... .......... .......... .......... .......... 96%  128M 0s
    ## 129100K .......... .......... .......... .......... .......... 96% 36.7M 0s
    ## 129150K .......... .......... .......... .......... .......... 96% 34.1M 0s
    ## 129200K .......... .......... .......... .......... .......... 96% 91.5M 0s
    ## 129250K .......... .......... .......... .......... .......... 96% 8.86M 0s
    ## 129300K .......... .......... .......... .......... .......... 96%  134M 0s
    ## 129350K .......... .......... .......... .......... .......... 96%  114M 0s
    ## 129400K .......... .......... .......... .......... .......... 96%  137M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  112M 0s
    ## 129500K .......... .......... .......... .......... .......... 96%  124M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 24.2M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 89.5M 0s
    ## 129650K .......... .......... .......... .......... .......... 96% 34.3M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 92.9M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 88.2M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 43.8M 0s
    ## 129850K .......... .......... .......... .......... .......... 96%  139M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 39.3M 0s
    ## 129950K .......... .......... .......... .......... .......... 96%  113M 0s
    ## 130000K .......... .......... .......... .......... .......... 97% 50.4M 0s
    ## 130050K .......... .......... .......... .......... .......... 97% 92.9M 0s
    ## 130100K .......... .......... .......... .......... .......... 97%  105M 0s
    ## 130150K .......... .......... .......... .......... .......... 97%  120M 0s
    ## 130200K .......... .......... .......... .......... .......... 97%  121M 0s
    ## 130250K .......... .......... .......... .......... .......... 97%  141M 0s
    ## 130300K .......... .......... .......... .......... .......... 97%  140M 0s
    ## 130350K .......... .......... .......... .......... .......... 97% 25.4M 0s
    ## 130400K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 130450K .......... .......... .......... .......... .......... 97% 54.0M 0s
    ## 130500K .......... .......... .......... .......... .......... 97%  129M 0s
    ## 130550K .......... .......... .......... .......... .......... 97% 77.7M 0s
    ## 130600K .......... .......... .......... .......... .......... 97%  106M 0s
    ## 130650K .......... .......... .......... .......... .......... 97%  105M 0s
    ## 130700K .......... .......... .......... .......... .......... 97%  117M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 82.5M 0s
    ## 130800K .......... .......... .......... .......... .......... 97%  124M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 87.6M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 56.9M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 89.0M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 83.3M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 80.5M 0s
    ## 131100K .......... .......... .......... .......... .......... 97%  105M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 58.7M 0s
    ## 131200K .......... .......... .......... .......... .......... 97%  116M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 98.0M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 83.8M 0s
    ## 131350K .......... .......... .......... .......... .......... 98% 69.0M 0s
    ## 131400K .......... .......... .......... .......... .......... 98%  129M 0s
    ## 131450K .......... .......... .......... .......... .......... 98% 95.3M 0s
    ## 131500K .......... .......... .......... .......... .......... 98% 72.3M 0s
    ## 131550K .......... .......... .......... .......... .......... 98% 70.6M 0s
    ## 131600K .......... .......... .......... .......... .......... 98% 85.2M 0s
    ## 131650K .......... .......... .......... .......... .......... 98% 74.2M 0s
    ## 131700K .......... .......... .......... .......... .......... 98% 93.3M 0s
    ## 131750K .......... .......... .......... .......... .......... 98% 86.4M 0s
    ## 131800K .......... .......... .......... .......... .......... 98%  105M 0s
    ## 131850K .......... .......... .......... .......... .......... 98% 70.7M 0s
    ## 131900K .......... .......... .......... .......... .......... 98%  104M 0s
    ## 131950K .......... .......... .......... .......... .......... 98% 73.9M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 75.5M 0s
    ## 132050K .......... .......... .......... .......... .......... 98% 94.7M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 92.2M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 92.2M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 86.9M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 83.4M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 96.7M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 13.8M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 74.7M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 74.1M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 80.6M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 94.1M 0s
    ## 132600K .......... .......... .......... .......... .......... 98%  117M 0s
    ## 132650K .......... .......... .......... .......... .......... 98%  152M 0s
    ## 132700K .......... .......... .......... .......... .......... 99% 75.6M 0s
    ## 132750K .......... .......... .......... .......... .......... 99% 44.4M 0s
    ## 132800K .......... .......... .......... .......... .......... 99% 76.5M 0s
    ## 132850K .......... .......... .......... .......... .......... 99%  108M 0s
    ## 132900K .......... .......... .......... .......... .......... 99% 50.4M 0s
    ## 132950K .......... .......... .......... .......... .......... 99% 98.9M 0s
    ## 133000K .......... .......... .......... .......... .......... 99%  121M 0s
    ## 133050K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 133100K .......... .......... .......... .......... .......... 99% 32.1M 0s
    ## 133150K .......... .......... .......... .......... .......... 99%  135M 0s
    ## 133200K .......... .......... .......... .......... .......... 99% 42.9M 0s
    ## 133250K .......... .......... .......... .......... .......... 99%  139M 0s
    ## 133300K .......... .......... .......... .......... .......... 99%  112M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  141M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 72.8M 0s
    ## 133450K .......... .......... .......... .......... .......... 99%  126M 0s
    ## 133500K .......... .......... .......... .......... .......... 99%  126M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 50.6M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 80.4M 0s
    ## 133650K .......... .......... .......... .......... .......... 99%  110M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 46.7M 0s
    ## 133750K .......... .......... .......... .......... .......... 99%  121M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 78.8M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 52.1M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 88.8M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 85.5M 0s
    ## 134000K .......... .......... .......... .......... .......... 99%  115M 0s
    ## 134050K .......... .....                                      100%  122M=3.6s
    ## 
    ## 2021-12-24 18:52:52 (36.1 MB/s) - ‘silva_nr99_v138.1_train_set.fa.gz.1’ saved [137283333/137283333]

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
```

##Assign taxonomy

``` r
fastaRef <-"/home/rstudio/silva_nr99_v138.1_train_set.fa.gz"
taxTab<-assignTaxonomy(seqtabNoC, refFasta=fastaRef, multithread=TRUE)
unname(head(taxTab))
```

    ##      [,1]       [,2]           [,3]          [,4]            [,5]            
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      [,6]         
    ## [1,] NA           
    ## [2,] NA           
    ## [3,] NA           
    ## [4,] NA           
    ## [5,] "Bacteroides"
    ## [6,] NA

##Construct phylogenetic tree

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

``` r
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

##Combine data into a phyloseq object

``` r
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy
all(rownames(seqtabAll) %in% samdf$SampleID) # TRUE
```

    ## [1] TRUE

``` r
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]
```

``` r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]

#Using phyloseq ##Loading the data

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

##Filtering ###Taxonomic Filtering

``` r
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

``` r
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

###Prevalence Filtering

``` r
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

##Agglomerate taxa

``` r
# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

##Abundance value transformation

``` r
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

``` r
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-41-1.png)<!-- --> ##Subset
by taxonomy

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->
###Preprocessing

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTTGGCAAGTTGGTAGCGTAATTCCGGGGCTCAACCTCGGCGCTACTACCAAAACTGCTGGACTTGAGTGCAGGAGGGGTGAATGGAATTCCTAGTGTAGCGGTGGAATGCGTAGATATTAGGAAGAACACCAGCGGCGAAGGCGATTCACTGGACTGTAACTGACACTGAGGCACGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->
#Different Ordination Projections

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTTGATAAGTCTGAAGTGAAAGGCCAAGGCTTAACCATGGAACTGCTTTGGAAACTATGAGGCTAGAGTGCTGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACAGAAACTGACGTTGAGGCTCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-53-1.png)<!-- --> ##Why are
the ordination plots so far from square? ###PCA on ranks

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

``` r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->
###Canonical correspondence

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGTGGCATGGCAAGCCAGAAGTGAAAACCCGGGGCTTAACCCCGCGGATTGCTTTTGGAACTGTCAGGCTGGAGTGCAGGAGGGGCAGGCGGAATTCCTGGTGTAGCGGTGAAATGCGTAGATATCAGGAGGAACACCGGTGGCGAAGGCGGCCTGCTGGACTGTAACTGACACTGAGGCTCGAAAGCGTGGGGA
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-EcoG2_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->
