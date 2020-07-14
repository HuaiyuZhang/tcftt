
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tcftt

<!-- badges: start -->

<!-- badges: end -->

The classical two-sample t-test only fit for the normal data. The tcfu()
and tt() tests implemented in this package are suitable for testing the
equality of two-sample means for the populations having unequal
variances. When the populations are not normally distributed, these
tests can provide more power than a large-sample t-test using normal
approximation, especially when the sample sizes are moderate. The tcfu()
uses the Cornish-Fisher expansion to achieve a better approximation to
the true percentiles. The tt() transforms the Welch’s t-statistic so
that the sampling distribution become more symmetric. More technical
details please refer to Zhang (2019) <http://hdl.handle.net/2097/40235>.

## Installation

You can install the released version of tcftt from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("tcftt")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tcftt)
x1 <- rnorm(20, 1, 3)
x2 <- rnorm(21, 2, 3)
tcfu(x1, x2, alternative = 'two.sided')
#> $stat
#> [1] -2.203019
#> 
#> $cutoff
#> [1] -1.974012  2.070452
#> 
#> $pvalue
#> [1] 0.0337871
#> 
#> $reject
#> [1] TRUE
tt(x1, x2, alternative = 'less')
#> $stat
#> [1] -2.262917
#> 
#> $cutoff
#> [1] -1.644854
#> 
#> $pvalue
#> [1] 0.9881796
#> 
#> $reject
#> [1] TRUE
```

## Main functions

The function `tcfu()` implements the Cornish-Fisher based two-sample
test (TCFU) and `tt()` implements the transformation based two-sample
test (TT).

The function `edgeworth()` provides the Edgeworth expansion for
cumulative density function the Welch’s t-statistic, and
`cornish_fisher()` provides the Cornish-Fisher expansion for the
percentiles.

The functions `adjust_power()` and `pauc()` provide power adjustment
methods for simulation studies.
