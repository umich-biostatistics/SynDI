
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `SynDI`

# Synthetic Data Integration

[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/umich-biostatistics/SynDI)
[![](https://img.shields.io/github/languages/code-size/umich-biostatistics/SynDI.svg)](https://github.com/umich-biostatistics/SynDI)

## Overview

Regression inference for multiple populations by integrating
summary-level data using stacked imputations. Gu, T., Taylor, J.M.G. and
Mukherjee, B. (2021) Regression inference for multiple populations by
integrating summary-level data using stacked imputations
\<arXiv:2106.06835>.

## Installation

If the devtools package is not yet installed, install it first:

``` r
install.packages('devtools')
```

``` r
# install the package from Github:
devtools::install_github('umich-biostatistics/SynDI', build_vignettes = TRUE) 
```

Once installed, load the package:

``` r
library(SynDI)
```

## Example Usage

For examples, see the package vignettes:

``` r
vignette("SynDI-example-binary")
vignette("SynDI-example-continuous")
```

### Current Suggested Citation

Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression inference for
multiple populations by integrating summary-level data using stacked
imputations arXiv:2106.06835.
