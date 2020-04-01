Distributions of Possible Effects
================

### Why should I use it?

The threat of endogeneity is ubiquitous within applied empirical
research. Regressor-error dependencies are a common inferential concern
not in the least because they may arise from any combination of omitted
variable, systematic measurement errors, self selection, systematic
missing data, reciprocal causation, or interference between units.
Conventional statistical methods do not reflect any source of
uncertainty other than random error and so do not express these
additional doubts. The package impliments a “near Bayesian” method of
sensitivity analysis which samples uniformly from the set of valid
control functions to build a distribution of possible causal effects as
well as graphical tools for assessing the sensitivity of one’s results
to the threat of hidden biases.

You can find a draft of the working paper introducing the approach

### How do I get it?

Until the package is released on CRAN you can install the developmental
version of the package with the following lines of code:

``` r
devtools::install_github("christophercschwarz/DOPE")
```

    ## Skipping install of 'DOPE' from a github remote, the SHA1 (a4448997) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(DOPE)
```

### How do I use it?

The DOPE algorithm is built around linear models
