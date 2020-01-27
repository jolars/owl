
<!-- README.md is generated from README.Rmd. Please edit that file -->

# owl <img src='man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jolars/owl.svg?branch=master)](https://travis-ci.org/jolars/owl)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/jolars/owl?branch=master&svg=true)](https://ci.appveyor.com/project/jolars/owl)
[![Coverage
status](https://codecov.io/gh/jolars/owl/branch/master/graph/badge.svg)](https://codecov.io/github/jolars/owl?branch=master)
<!-- badges: end -->

Efficient implementations for Sorted L-One Penalized Estimation (SLOPE):
generalized linear models regularized with the sorted L1-norm or,
equivalently, ordered weighted L1-norm (OWL). There is support for
ordinary least-squares regression, binomial regression, multinomial
regression, and poisson regression, as well as both dense and sparse
predictor matrices. In addition, the package features predictor
screening rules that enable efficient solutions to high-dimensional
problems.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("jolars/owl")
```

## License

owl is open source software, licensed under
[GPL-3](https://github.com/jolars/owl/blob/master/LICENSE).

## Versioning

owl uses [semantic versioning](http://semver.org).

## Code of conduct

Please note that **owl** is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree
to abide by its terms.
