
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tcftt

<!-- badges: start -->

<!-- badges: end -->

The classical two-sample t-test only fits for the normal data. The
tcfu() and tt() tests implemented in this package are suitable for
testing the equality of two-sample means for the populations having
unequal variances. When the populations are not normally distributed,
these tests can provide more power than a large-sample t-test using
normal approximation, especially when the sample sizes are moderate. The
tcfu() uses the Cornish-Fisher expansion to achieve a better
approximation to the true percentiles. The tt() transforms the Welch’s
t-statistic so that the sampling distribution become more symmetric.
More technical details please refer to [Zhang (2019)](http://hdl.handle.net/2097/40235) and [Zhang and Wang (2021)](https://www.tandfonline.com/doi/abs/10.1080/10485252.2021.1982938).

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
#> [1] -1.044103
#> 
#> $cutoff
#> [1] -1.970350  2.073316
#> 
#> $pvalue
#> [1] 0.3019628
#> 
#> $reject
#> [1] FALSE
tt(x1, x2, alternative = 'less')
#> $stat
#> [1] -1.063013
#> 
#> $cutoff
#> [1] -1.644854
#> 
#> $pvalue
#> [1] 0.8561119
#> 
#> $reject
#> [1] FALSE
```

## Main functions

The function `tcfu()` implements the Cornish-Fisher based two-sample
test (TCFU) and `tt()` implements the transformation based two-sample
test (TT).

The function `t_edgeworth()` provides the Edgeworth expansion of the
cumulative density function for the Welch’s t-statistic, and
`t_cornish_fisher()` provides the Cornish-Fisher expansion for its
percentiles.

The functions `adjust_power()` and `pauc()` provide power adjustment
methods for simulation studies.
