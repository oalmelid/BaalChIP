We thank the [Breast Cancer Research Foundation](https://give.bcrfcure.org/checkout/donation?eid=31404&_ga=1.141086386.1225939115.1465824174) for financial support of this work.

# Installation

To install the development version of BaalChIP R package from GitHub type:

```r
library(devtools)
install_github("inesdesantiago/BaalChIP")
```

# BaalChIP R package

BaalChIP (Baysian Anaysis of Allelic imbalances from ChIP-seq data) was built for the identification of allele-specific binding (ASB) events from ChIP-seq data obtained from cancer cell lines.

BaalChIP tests the differential read counts of the alleles at each heterozygous variant using the quantitative information of ChIP-seq read counts at the reference and alternative alleles and accomodating the information about the allele presence and other sources of ChIP-seq mapping biases.

## Example run
This document offers a quick example of how to use BaalChIP to identify ASB events with correction for relative allele frequency.

The example dataset contains ChIP-seq data obtained for two cell lines: A cancer cell-line (MCF7) and a normal cell line (GM12891). For each cell line, ChIP-seq data exists for four transcription factors and two biological replicates for each of the transcription factors.
The metadata and all files necessary for this example are available in the extra subdirectory of the BaalChIP package directory; you can make this your working directory by entering:

```r
library(BaalChIP)
setwd(system.file("test",package="BaalChIP"))
```
Note that the example data in this vignette does not reveal real biology and was build only for demonstration purposes. 

The first step is to contruct a BaalChIP object:

```r
samplesheet <- "exampleChIP.tsv"
hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
```
Given a new BaalChIP object, allele-specific binding events can be identified as follows:

```r
#first load some data
data(blacklist_hg19)
data(pickrell2011cov1_hg19)
data(UniqueMappability50bp_hg19)

#run example
res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
res <- QCfilter(res,
                RegionsToFilter=list("blacklist"=blacklist_hg19,
                                     "highcoverage"=pickrell2011cov1_hg19),
                RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
res <- mergePerGroup(res)
res <- filter1allele(res)
res <- getASB(res, Iter=5000, conf_level=0.95)
```

## Example run in a 2-step wrapper script

If you trust the package defaults, the first four steps can be replaced by a wrapper function, making BaalChIP workflow possible to run in a 3-step script:

```r
res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
res <- BaalChIP.run(res)
```

The package vignette describes these steps in more detail.
