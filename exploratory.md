aortopathy
================
John Chamberlin
10/19/2022

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.1.2

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## Warning: package 'tidyr' was built under R version 4.1.2

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 4.1.2

    ## ── Attaching packages
    ## ───────────────────────────────────────
    ## tidyverse 1.3.2 ──

    ## ✓ ggplot2 3.3.6     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.2

    ## Warning: package 'ggplot2' was built under R version 4.1.2

    ## Warning: package 'readr' was built under R version 4.1.2

    ## Warning: package 'forcats' was built under R version 4.1.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(tximport)
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.1.3

    ## Loading required package: stats4
    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min
    ## 
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 4.1.2

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.1.2

    ## Loading required package: SummarizedExperiment
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
# import the metadata
m1 = readxl::read_xlsx("data/MetadataSamples_QuinlanLab_5.27.xlsx")
m2 = readxl::read_xlsx("data/Glotzbach Sample Match List 05312022.xlsx")

# we care about the GLO ID, ID and the phenotype
m1 = m1[complete.cases(m1),] # drop empty row
m2 = m2 %>% select_if(function(x) !(all(is.na(x))))
m2 = m2 %>% filter(!is.na(Organism))

colnames(m1) = make.names(colnames(m1),unique=TRUE)
colnames(m2) = make.names(colnames(m2),unique=TRUE)

print((m1 %>% left_join(m2, by = c("Glotzbach.Lab.ID" = "Sample.Name"))))
```

    ## # A tibble: 20 × 20
    ##    Phenotype         Glotzbach.Lab.ID HCI.EXP. Storage Conc...ng.uL. Sample.Type
    ##    <chr>             <chr>            <chr>    <chr>   <chr>         <chr>      
    ##  1 Bicuspid valve +… GLO-033          18845X1  RNAlat… 29.2          Total RNA …
    ##  2 Bicuspid valve +… GLO-038          18845X2  RNAlat… 22.4          Total RNA …
    ##  3 Bicuspid valve +… GLO-039          18845X3  RNAlat… 40            Total RNA …
    ##  4 Bicuspid valve +… GLO-046          18845X4  RNAlat… 31.3          Total RNA …
    ##  5 Tricuspid valve … GLO-045          18845X5  RNAlat… 8.6999999999… Total RNA …
    ##  6 Tricuspid valve … GLO-037          18845X6  RNAlat… 44.7          Total RNA …
    ##  7 Donor: tricuspid… GLO-041          18845X7  RNAlat… 40.4          Total RNA …
    ##  8 Donor: tricuspid… GLO-041          18845X7  RNAlat… 52.7          Total RNA …
    ##  9 Donor: tricuspid… GLO-044          18845X8  RNAlat… 28.8          Total RNA …
    ## 10 Donor: tricuspid… GLO-044          18845X8  RNAlat… 32.1          Total RNA …
    ## 11 Bicuspid valve +… GLO-031          18794X1  Fresh … <NA>          <NA>       
    ## 12 Bicuspid valve +… GLO-008          18794X2  Fresh … <NA>          <NA>       
    ## 13 Bicuspid valve +… GLO-033          18794X3  Fresh … 29.2          Total RNA …
    ## 14 Bicuspid valve +… GLO-013          18794X6  Fresh … <NA>          <NA>       
    ## 15 Bicuspid valve +… GLO-024          18794X7  Fresh … <NA>          <NA>       
    ## 16 Bicuspid valve +… GLO-028          18794X10 Fresh … <NA>          <NA>       
    ## 17 Bicuspid valve +… GLO-032          18794X9  Fresh … <NA>          <NA>       
    ## 18 Tricuspid valve … GLO-027          18794X5  Fresh … <NA>          <NA>       
    ## 19 Tricuspid valve … GLO-004          18794X8  Fresh … <NA>          <NA>       
    ## 20 Donor: tricuspid… GLO-003          18794X4  Fresh … <NA>          <NA>       
    ## # … with 14 more variables: Organism <chr>, QC.Conc...ng.uL. <dbl>,
    ## #   QC.RIN <dbl>, Library.Protocol <chr>, Prepped.by.Core. <chr>,
    ## #   Index.Tag.A <chr>, Index.Tag.Sequence.A <chr>, Index.Tag.B <chr>,
    ## #   Index.Tag.Sequence.B <chr>, ID <chr>, Lib.QC.Conc. <dbl>, QC.Status <chr>,
    ## #   Seq.Lib.Prep.Status <chr>, Core.to.prep.lib. <chr>

``` r
ifiles = list.files("data/rsem/",
                    pattern = "*.isoforms.results")
gfiles = list.files("data/rsem/",
                    pattern = "*.genes.results")

files = paste0("data/rsem/",ifiles)
names(files) = gsub("[:_:].*$","",ifiles)

tx2gene = data.table::fread("data/rsem/18890X1_210521_A00421_0326_BH3HT3DSX2_S1_L001.isoforms.results") %>% select(transcript_id, gene_id)

txi.rsem = tximport(files,type = "rsem", 
                    tx2gene = tx2gene)
```

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
sampleTable = data.frame("batch"=gsub("X.*$","",colnames(txi.rsem$counts)))
rownames(sampleTable) = colnames(txi.rsem$counts)

dds <- DESeq2::DESeqDataSetFromTximport(txi.rsem, sampleTable, ~batch)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
vsd <- vst(dds, blind=FALSE)
```

    ## using 'avgTxLength' from assays(dds), correcting for library size

``` r
plotPCA(vsd, intgroup = "batch")
```

![](exploratory_files/figure-gfm/pca-1.png)<!-- -->

``` r
# toy example of differential expression between batches, not phenotype
dds = DESeq(dds)
```

    ## estimating size factors

    ## using 'avgTxLength' from assays(dds), correcting for library size

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 458 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
res = results(dds)
res
```

    ## log2 fold change (MLE): batch 19192 vs 18890 
    ## Wald test p-value: batch 19192 vs 18890 
    ## DataFrame with 60664 rows and 6 columns
    ##                    baseMean log2FoldChange     lfcSE       stat      pvalue
    ##                   <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ## ENSG00000000003    134.1600      0.5161879 0.1792438   2.879809 3.97916e-03
    ## ENSG00000000005     11.4353      1.3906709 1.6798169   0.827870 4.07744e-01
    ## ENSG00000000419    395.7892     -1.0085819 0.0944267 -10.681114 1.24765e-26
    ## ENSG00000000457    210.8617      0.0703597 0.1004840   0.700207 4.83798e-01
    ## ENSG00000000460     76.8214     -0.3251698 0.1935670  -1.679882 9.29802e-02
    ## ...                     ...            ...       ...        ...         ...
    ## ENSG00000288721  20.3835160       1.030731  0.335484   3.072372 2.12365e-03
    ## ENSG00000288722 300.4529643      -0.490747  0.289480  -1.695272 9.00238e-02
    ## ENSG00000288723   0.0531436       1.059805  3.123483   0.339302 7.34382e-01
    ## ENSG00000288724   0.0000000             NA        NA         NA          NA
    ## ENSG00000288725  91.8659734      -5.431470  1.318161  -4.120490 3.78068e-05
    ##                        padj
    ##                   <numeric>
    ## ENSG00000000003 1.35378e-02
    ## ENSG00000000005 5.72933e-01
    ## ENSG00000000419 6.45287e-25
    ## ENSG00000000457 6.42981e-01
    ## ENSG00000000460 1.91408e-01
    ## ...                     ...
    ## ENSG00000288721 0.007775382
    ## ENSG00000288722 0.186494804
    ## ENSG00000288723          NA
    ## ENSG00000288724          NA
    ## ENSG00000288725 0.000204751

``` r
plotMA(res, ylim = c(-3,3))
```

![](exploratory_files/figure-gfm/difexp-1.png)<!-- -->
