# SNPClust

<!--
[![Build Status](https://travis-ci.org/yihui/knitr.svg)](https://travis-ci.org/yihui/knitr)
[![Coverage Status](https://coveralls.io/repos/yihui/knitr/badge.svg?branch=master&service=github)](https://coveralls.io/github/yihui/knitr?branch=master)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/knitr)](http://cran.rstudio.com/package=knitr)
-->

The R package **snpclust** performs unsupervised feature selection and summarization based on Principal Component Analysis, Gaussian Mixture Models, and Markov Chain Monte Carlo.

## Installation

SNPClust is an R package that requires the PLINK and SHAPEIT softwares.
You can use the **devtools** R package to install the development version from Github:

```r
devtools::install_github('ThomasChln/snpclust', build_vignettes = TRUE)
```

## Motivation

SNPClust is developed for the reclassification of Systemic Autoimmune Diseases based on genetic markers instead of clinical criteria for the European-funded project [PRECISESADS](http://precisesads.eu).

## Usage

The snpclust function is the main function and can be called on a Genomic Data Structure file path.
Results can be displayed with several plots. A demo is available in the *SNPClust applied to Europeans* HTML vignette.

## License

This package is free and open source software, licensed under GPL-3.
