# SNPClust

<!--
[![Build Status](https://travis-ci.org/yihui/knitr.svg)](https://travis-ci.org/yihui/knitr)
[![Coverage Status](https://coveralls.io/repos/yihui/knitr/badge.svg?branch=master&service=github)](https://coveralls.io/github/yihui/knitr?branch=master)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/knitr)](http://cran.rstudio.com/package=knitr)
-->

The R package **snpclust** performs unsupervised feature selection and summarization based on Principal Component Analysis, Gaussian Mixture Models, and Markov Chain Monte Carlo.

## Installation

SNPClust is an [R](https://cran.r-project.org/) package that requires the softwares [PLINK](http://zzz.bwh.harvard.edu/plink/download.shtml) and [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download).
Once they are installed, you can open R and use the **devtools** R package to install the development version from Github:
```r
if (!require(devtools)) install.packages('devtools')
devtools::install_github('ThomasChln/snpclust', build_vignettes = TRUE)
```

## Motivation

SNPClust is developed for the reclassification of Systemic Autoimmune Diseases (SADs) based on genetic markers instead of clinical criteria for the European-funded project [PRECISESADS](http://precisesads.eu).

The results of the application to SADs were presented in the [Genome-wide unsupervised clustering](http://f1000research.com/posters/1098306) poster at the [Basel Computational Biology Conference](http://www.bc2.ch) in June 2015.

## Usage

The snpclust function is the main function and can be called on a snpgds file path from the [R package SNPRelate](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html).
Results can be displayed with several plots. A demo is available in the *SNPClust applied to Europeans* HTML vignette.

## License

This package is free and open source software, licensed under GPL-3.
