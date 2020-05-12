# SNPClust

[![Build Status](https://travis-ci.org/thomaschln/snpclust.svg)](https://travis-ci.org/thomaschln/snpclust)
[![Coverage status](https://codecov.io/gh/thomaschln/snpclust/branch/master/graph/badge.svg)](https://codecov.io/github/thomaschln/snpclust)

The R package **snpclust** performs unsupervised feature selection and summarization by selecting SNPs based on principal component analysis (PCA) contributions, and estimating haplotypes of nearby correlated SNPs using the SHAPEIT software.

## Installation

### Devtools

SNPClust is an [R](https://cran.r-project.org/) package that requires the softwares [PLINK](http://zzz.bwh.harvard.edu/plink/download.shtml) and [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download).
Once they are installed, you can open R and use the **devtools** R package to install the development version from Github:
```r
if (!require(devtools)) install.packages('devtools')
devtools::install_github('ThomasChln/snpclust', build_vignettes = TRUE)
```

### Docker

A Docker image with PLINK, SHAPEIT, and SNPClust installed is available.

```
docker pull thomaschln/snpclust
```

## Motivation

SNPClust was developed to reclassify systemic autoimmune diseases (SADs) based on genetic markers instead of clinical criteria for the European-funded project PRECISESADS.

It was applied to a genome wide dataset of 379,190 SNPs from 4,212 systemic lupus erythematosus (SLE) patients and 1,221 healthy controls and results were published in PLOS ONE: [Single Nucleotide Polymorphism Clustering in Systemic Autoimmune Diseases](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160270)

It was also used to reproduce the Human Genome Diversity Panel PCA of 300,000 SNPs of 1,000 samples [https://f1000research.com/articles/6-278/v1].

## Usage

The snpclust function is the main function and can be called on a snpgds file path from the [R package SNPRelate](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html).
Results can be displayed with several plots. A demo is available in the *SNPClust applied to Europeans* HTML vignette.

## License

This package is free and open source software, licensed under GPL-3.
