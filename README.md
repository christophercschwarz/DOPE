Distributions of Possible Effects
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

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
control functions under ignorance to build a distribution of possible
causal effects. It also provides graphical tools for assessing the
sensitivity of one’s results to the threat of hidden biases.

### What does it do?

Regardless of the source a regressor-error dependency is a
regressor-error dependency and can be dealt with in the same manner;
namely, via control functions. Unlike a number of existing forms of
sensitivity analysis which treat omitted variables as one of many
sources of endogeneity, the control function approach treats all
endogeneity as an omitted variables problem.

A simple example aids intuition. Suppose that we are interested in the
effect of some treatment(s)
![\\mathbf{X}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BX%7D
"\\mathbf{X}") on outcome
![\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BY%7D
"\\mathbf{Y}"). Assuming constant effects, linearity, and including the
intercept in the design matrix we can use the common linear
specification

<center>

  
![\\mathbf{Y} = \\mathbf{X\\beta} +
\\mathbf{\\epsilon^\*}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BY%7D%20%3D%20%5Cmathbf%7BX%5Cbeta%7D%20%2B%20%5Cmathbf%7B%5Cepsilon%5E%2A%7D
"\\mathbf{Y} = \\mathbf{X\\beta} + \\mathbf{\\epsilon^*}")  

</center>

where
![\\mathbf{\\epsilon^\*}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7B%5Cepsilon%5E%2A%7D
"\\mathbf{\\epsilon^*}") is a latent random variable and
![\\beta](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta
"\\beta") is the ATE. Regardless of the source, regressor-error
dependencies manifest such that
![E(\\mathbf{\\epsilon^\*}|\\mathbf{X})=0](http://chart.apis.google.com/chart?cht=tx&chl=E%28%5Cmathbf%7B%5Cepsilon%5E%2A%7D%7C%5Cmathbf%7BX%7D%29%3D0
"E(\\mathbf{\\epsilon^*}|\\mathbf{X})=0") and so

<center>

![\\begin{align} \\hat{\\mathbf{\\beta}} &=
(\\mathbf{X}^T\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{Y} \\\\ &=
\\mathbf{\\beta} +
(\\mathbf{X}^T\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{\\epsilon^\*} \\neq
\\beta
\\end{align}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbegin%7Balign%7D%20%20%5Chat%7B%5Cmathbf%7B%5Cbeta%7D%7D%20%26%3D%20%28%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BX%7D%29%5E%7B-1%7D%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BY%7D%20%5C%5C%20%20%26%3D%20%5Cmathbf%7B%5Cbeta%7D%20%2B%20%28%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BX%7D%29%5E%7B-1%7D%5Cmathbf%7BX%7D%5ET%5Cmathbf%7B%5Cepsilon%5E%2A%7D%20%5Cneq%20%5Cbeta%20%5Cend%7Balign%7D
"\\begin{align}  \\hat{\\mathbf{\\beta}} &= (\\mathbf{X}^T\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{Y} \\\\  &= \\mathbf{\\beta} + (\\mathbf{X}^T\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{\\epsilon^*} \\neq \\beta \\end{align}")

</center>

The basic idea of a control function is to construct, usually using
instrumental variables, some
![\\mathbf{V}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BV%7D
"\\mathbf{V}") for which
![E(\\mathbf{\\epsilon^\*}|\\mathbf{X},\\mathbf{V})=0](http://chart.apis.google.com/chart?cht=tx&chl=E%28%5Cmathbf%7B%5Cepsilon%5E%2A%7D%7C%5Cmathbf%7BX%7D%2C%5Cmathbf%7BV%7D%29%3D0
"E(\\mathbf{\\epsilon^*}|\\mathbf{X},\\mathbf{V})=0") and include that
term as a regressor to purge the dependency. We are usually in the
position where no such control function can be constructed. Regardless,
we know that with the correct control function the regression of
![\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BY%7D
"\\mathbf{Y}") on
![\\mathbf{Z}=(\\mathbf{X},\\mathbf{V})](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BZ%7D%3D%28%5Cmathbf%7BX%7D%2C%5Cmathbf%7BV%7D%29
"\\mathbf{Z}=(\\mathbf{X},\\mathbf{V})") yields coefficients

<center>

  
![\\hat{\\mathbf{\\beta}} =
(\\mathbf{Z}^T\\mathbf{Z})^{-1}\\mathbf{Z}^T\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cmathbf%7B%5Cbeta%7D%7D%20%3D%20%28%5Cmathbf%7BZ%7D%5ET%5Cmathbf%7BZ%7D%29%5E%7B-1%7D%5Cmathbf%7BZ%7D%5ET%5Cmathbf%7BY%7D
"\\hat{\\mathbf{\\beta}} = (\\mathbf{Z}^T\\mathbf{Z})^{-1}\\mathbf{Z}^T\\mathbf{Y}")  

</center>

which contains coefficients carrying a causal interpretation.

So what we will do is sample uniformly from the set of valid control
functions to build out the distribution of causal effects given our
ignorance. An algorithm to generate one such randomly drawn control
function is the following:

<center>

<a href="https://imgbb.com/"><img src="https://i.ibb.co/Wng5xrm/algorithm.png" alt="algorithm" border="0"></a>

</center>

This algorithm effectively repeatedly uses the Cholesky decomposition to
build up a valid augmented covariance matrix from which quantities of
interest can be calculated.

To apply this to regression, note that we may rearrange the augmented
covariance matrix as

<center>

![\\Sigma^\* = \\begin{bmatrix}
\\mathbf{K}\_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}} &
\\mathbf{K}\_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}\\\\
\\mathbf{K}\_{\\mathbf{Y}\\tilde{\\mathbf{Z}}} & \\mathbf{K\_{YY}}
\\end{bmatrix}](http://chart.apis.google.com/chart?cht=tx&chl=%5CSigma%5E%2A%20%3D%20%5Cbegin%7Bbmatrix%7D%20%20%5Cmathbf%7BK%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%7D%20%26%20%5Cmathbf%7BK%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Cmathbf%7BY%7D%7D%5C%5C%20%20%5Cmathbf%7BK%7D_%7B%5Cmathbf%7BY%7D%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%7D%20%26%20%5Cmathbf%7BK_%7BYY%7D%7D%20%5Cend%7Bbmatrix%7D
"\\Sigma^* = \\begin{bmatrix}  \\mathbf{K}_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}} & \\mathbf{K}_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}\\\\  \\mathbf{K}_{\\mathbf{Y}\\tilde{\\mathbf{Z}}} & \\mathbf{K_{YY}} \\end{bmatrix}")

</center>

with ![\\tilde{\\mathbf{Z}} =
(\\mathbf{X},\\tilde{\\mathbf{V}})](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%20%3D%20%28%5Cmathbf%7BX%7D%2C%5Ctilde%7B%5Cmathbf%7BV%7D%7D%29
"\\tilde{\\mathbf{Z}} = (\\mathbf{X},\\tilde{\\mathbf{V}})") where
![\\tilde{\\mathbf{V}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7B%5Cmathbf%7BV%7D%7D
"\\tilde{\\mathbf{V}}") is the control function defined by one run of
Algorithm 1. The coefficients from the regression of
![\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BY%7D
"\\mathbf{Y}") on
![\\mathbf{X}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BX%7D
"\\mathbf{X}") and
![\\tilde{\\mathbf{V}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7B%5Cmathbf%7BV%7D%7D
"\\tilde{\\mathbf{V}}") are given by

<center>

  
![\\tilde{\\mathbf{\\beta}} =
\\mathbf{K}^{-1}\_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}}\\mathbf{K}\_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7B%5Cmathbf%7B%5Cbeta%7D%7D%20%3D%20%5Cmathbf%7BK%7D%5E%7B-1%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%7D%5Cmathbf%7BK%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Cmathbf%7BY%7D%7D
"\\tilde{\\mathbf{\\beta}} = \\mathbf{K}^{-1}_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}}\\mathbf{K}_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}")  

</center>

so one never actually have to generate a realization of
![\\tilde{V}](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7BV%7D
"\\tilde{V}") once one has the covariances by which it is defined. We
may similarly derive the R-squared with this information alone:

<center>

![\\tilde{R}^2 =
\\frac{\\mathbf{K}^{T}\_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}\\mathbf{K}^{-1}\_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}}\\mathbf{K}\_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}}{\\mathbf{K}\_{\\mathbf{Y}\\mathbf{Y}}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Ctilde%7BR%7D%5E2%20%3D%20%5Cfrac%7B%5Cmathbf%7BK%7D%5E%7BT%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Cmathbf%7BY%7D%7D%5Cmathbf%7BK%7D%5E%7B-1%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%7D%5Cmathbf%7BK%7D_%7B%5Ctilde%7B%5Cmathbf%7BZ%7D%7D%5Cmathbf%7BY%7D%7D%7D%7B%5Cmathbf%7BK%7D_%7B%5Cmathbf%7BY%7D%5Cmathbf%7BY%7D%7D%7D
"\\tilde{R}^2 = \\frac{\\mathbf{K}^{T}_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}\\mathbf{K}^{-1}_{\\tilde{\\mathbf{Z}}\\tilde{\\mathbf{Z}}}\\mathbf{K}_{\\tilde{\\mathbf{Z}}\\mathbf{Y}}}{\\mathbf{K}_{\\mathbf{Y}\\mathbf{Y}}}")

</center>

By drawing thousands of these control functions we get thousands of
coefficient vectors and R-squared values which represent the additional
uncertainty in our results to the threat of hidden biases.

### How do I get it?

Until the package is released on CRAN you can install the developmental
version of the package with the following lines of code:

``` r
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
devtools::install_github("christophercschwarz/DOPE",
                         dependencies=TRUE)
library(DOPE)
```

### How do I use it?

The DOPE algorithm is built upon linear regression targeting uncertainty
in the ATE. The approach can be extended to semi-parametric
distributional regression models, possibly including random effects, for
various estimands with relative ease utilizing whitening and
augmentation tricks to put the models in OLS form.

To illustrate a number of functions from the package, let’s look at a
simulated example. First, let’s generate a random correlation matrix,
generate some data, and run a linear model.

``` r
set.seed(8675309)
x_vars <- 5
n_obs <- 1000
corm <- RandomCormCPP(nvars = x_vars)
X_mat <- MASS::mvrnorm(n_obs, rep(0,x_vars), Sigma = corm)

betas <- 1:x_vars

y <- X_mat %*% betas + rnorm(n_obs, 0, 5)

dat <- data.frame(y,X_mat)
cov(dat)
##             y         V1         V2         V3         V4         V5
## y  121.467197  8.2838206 -3.3630478  6.0726558  8.1503611  8.5835385
## V1   8.283821  0.9339341 -0.6078563  0.5360767  0.7143011  0.7979034
## V2  -3.363048 -0.6078563  0.9102019 -0.6559167 -0.1451956 -0.3767998
## V3   6.072656  0.5360767 -0.6559167  1.0287591  0.3094510  0.5009242
## V4   8.150361  0.7143011 -0.1451956  0.3094510  1.0035173  0.5555419
## V5   8.583538  0.7979034 -0.3767998  0.5009242  0.5555419  0.9444915
```

The lower V1-V5 sub-matrix is a draw from the set of valid correlation
matrices by repeatedly applying the DOPE augmentation algorithm. Within
the context of regression, what we will do is take the observed
covariance matrix as given and draw an additional row and column from
the set of valid covariance matrices, calculate quantities of interest,
and repeat. To begin, let’s estimate the linear model:

``` r
mod <- lm(y ~ ., data=dat)
```

Since we have control over the process by which the data was generated
we know that the model is correct and that the coefficients are
unbiased/consistent estimates of the treatment effect of interest.
Suppose, however, that we did not. We can take draws from the set of
valid control functions to build out a distribution of possible effects
reflecting ignorance of the regressor-error dependency plaguing our
analysis. The number of draws can be set with `nsims` and the number of
cores for parallel computation with `n.cores`.

``` r
set.seed(6161918)
dope <- DOPE(mod, nsims = 5000, n.cores = parallel::detectCores())
```

``` r
tail(dplyr::tbl_df(dope))
## # A tibble: 6 x 8
##   Intercept     V1    V2     V3    V4      V5 ControlFunction R_Squared
##       <dbl>  <dbl> <dbl>  <dbl> <dbl>   <dbl>           <dbl>     <dbl>
## 1    0.385  28.1    7.96 -3.76  -2.10    7.36         -18.6       0.974
## 2    0.0843  0.830  1.57  2.90   4.16    5.18           0.243     0.788
## 3    0.0492  7.61   5.56  3.88  13.7    -1.62         -12.0       0.975
## 4    0.312  16.1    1.52 -0.900 -1.35    1.64          -9.09      0.999
## 5   -1.20   30.7   23.7  60.6   21.9  -157.           113.        0.846
## 6    0.0747  0.628  1.56  2.75   4.12    5.30          NA         0.788
```

The result is a dataframe of `nsims` + 1 observations with columns for
each of the estimated coefficients, the control function coefficient,
and model R-squared. The last observation simply re-states the results
from the naive model for use in plotting functions. The most basic plot
shows the simulated distribution of possible effects and gives a number
of useful summaries. Because these are `ggplot` objects, that can be
easily modified by adding additional layers.

``` r
plot_DOPE(dope,"V2") + ggtitle("Distribution of Possible Effects: V2")
```

<img src="README_figs/README-unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

In this example, based upon 5000 draws around 60% of the estimated
effects are greater than zero with a 95% credible interval given by the
lower and upper 95 bounds. The naive estimate is indicated in red. Since
the distribution is simulated setting a seed and using a large number of
draws is highly recommended, but these quantities usually become stable
with over 20,000 draws.

We can get a little bit more information from the plot by shading based
upon the R-squared from each regression and turning off the naive
result.

``` r
plot_DOPE(dope,"V2",shade=T,include_naive = F) + ggtitle("Distribution of Possible Effects: V2")
```

<img src="README_figs/README-unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

Lighter shades indicate lower R-squared than darker shades, given a
sense for how fits are distributed across the effect values. This is
important as the distribution of possible effects changes as a function
of the allowable R-squared. If we wanted to change the color away from
greyscale we can add another layer like so:

``` r

plot_DOPE(dope,"V2",shade=T,include_naive = F) + ggtitle("Distribution of Possible Effects: V2") + scale_fill_discrete()
```

<img src="README_figs/README-unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

If there is a particular color palette you are looking for you can use
it using something like the following.

``` r
library(wesanderson)

p <- plot_DOPE(dope,"V2",shade=T,include_naive = F) + ggtitle("Distribution of Possible Effects: V2") 

file <- tempfile()
st <- try(ggsave(file,device="pdf",p + scale_fill_manual(values = "black")),silent=T)
unlink(file)
n_fills <- as.numeric(regmatches(st,gregexpr("[[:digit:]]+",st))[[1]][1])

p + scale_fill_manual(values=colorRampPalette(wes_palette("Zissou1"))(n_fills))
```

<img src="README_figs/README-unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

We can get a sense for how the distribution of possible effects changes
with restrictions on the R-squared with a `sensitivity_plot`.

``` r
sensitivity_plot(dope,"V2")
```

<img src="README_figs/README-unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

The solid dots indicate how the certainty in the effect changes as you
say “the world is not that deterministic,” reducing the maximum
allowable R-squared. This reduces uncertainty relative to ignornace,
indicated by the dashed line, at the cost of restrictiveness, indicated
by the hollow points. In the other direction, one might say that “the
world is more deterministic than reflected in my model.” This is
represented by the crossed points and increases uncertainty until the
lower pessimistic bound on effect certainty. The proportion of draws
rejected by this lower thresholding is given by the filled diamonds.

In this particular example, there is no regressor-error dependency
actually biasing the result. The true R-squared is near the estimated
R-squared, and accordingly when we restrict attention to those estimates
close to this level we become certain that our estimates are
appropriately signed. In actual empirical practice, however, such
information is rarely known with any degree of confidence and so this
extra uncertainty should be reported.

### Notes on Development Path

There are a number of additional features that will be added to the
package over time. Of particular interest is adding additional
facilities for converting regression results into a format which may be
easily used with the package. Currently only linear models and GLMs are
supported, the latter requiring to be fit with the `DOPE_irls` function
which retains the IRLS weights and then subsequently “whitens” the data
to be fit with `lm`.

By staying firmly within a regression framework the package can be
extended to conduct sensitivity analysis on a wite array of
semiparametric distributional regression models as well as sampling
variability. For example, a bootstrap step could be implimented by
noting that ![n](http://chart.apis.google.com/chart?cht=tx&chl=n "n")
samples from ![n](http://chart.apis.google.com/chart?cht=tx&chl=n "n")
observations can be written as an ![n \\times
n](http://chart.apis.google.com/chart?cht=tx&chl=n%20%5Ctimes%20n
"n \\times n") matrix
![\\mathbf{Q}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BQ%7D
"\\mathbf{Q}"). For observed data
![\\mathbf{R}=(\\mathbf{Y},\\mathbf{X})](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BR%7D%3D%28%5Cmathbf%7BY%7D%2C%5Cmathbf%7BX%7D%29
"\\mathbf{R}=(\\mathbf{Y},\\mathbf{X})") a bootstrap replicate is simply
![\\mathbf{R}^\*=\\mathbf{QR}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BR%7D%5E%2A%3D%5Cmathbf%7BQR%7D
"\\mathbf{R}^*=\\mathbf{QR}"), the least squares estimate being

<center>

![\\beta^\* = (\\mathbf{X}^T \\mathbf{Q}^T \\mathbf{Q} \\mathbf{X})^{-1}
\\mathbf{X}^T \\mathbf{Q}^T
\\mathbf{QY}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta%5E%2A%20%3D%20%28%5Cmathbf%7BX%7D%5ET%20%5Cmathbf%7BQ%7D%5ET%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7BX%7D%29%5E%7B-1%7D%20%5Cmathbf%7BX%7D%5ET%20%5Cmathbf%7BQ%7D%5ET%20%5Cmathbf%7BQY%7D
"\\beta^* = (\\mathbf{X}^T \\mathbf{Q}^T \\mathbf{Q} \\mathbf{X})^{-1} \\mathbf{X}^T \\mathbf{Q}^T \\mathbf{QY}")

</center>

that is, a weighted least squares estimate. We can leverage the notion
of whitening from Aitken estimators like generalized least squares

<center>

  
![\\hat{\\beta} = (\\mathbf{X}^T \\Omega^{-1} \\mathbf{X})^T
\\mathbf{X}^T \\Omega^{-1}
\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cbeta%7D%20%3D%20%28%5Cmathbf%7BX%7D%5ET%20%5COmega%5E%7B-1%7D%20%5Cmathbf%7BX%7D%29%5ET%20%5Cmathbf%7BX%7D%5ET%20%5COmega%5E%7B-1%7D%20%5Cmathbf%7BY%7D
"\\hat{\\beta} = (\\mathbf{X}^T \\Omega^{-1} \\mathbf{X})^T \\mathbf{X}^T \\Omega^{-1} \\mathbf{Y}")  

</center>

and linearly transform the data where
![\\Omega](http://chart.apis.google.com/chart?cht=tx&chl=%5COmega
"\\Omega") is Cholesky decomposed into
![\\mathbf{CC}^T](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BCC%7D%5ET
"\\mathbf{CC}^T") and one regresses ![\\mathbf{Y}^\* = \\mathbf{C}^{-1}
\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BY%7D%5E%2A%20%3D%20%5Cmathbf%7BC%7D%5E%7B-1%7D%20%5Cmathbf%7BY%7D
"\\mathbf{Y}^* = \\mathbf{C}^{-1} \\mathbf{Y}") on ![\\mathbf{X}^\* =
\\mathbf{C}^{-1}
\\mathbf{X}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7BX%7D%5E%2A%20%3D%20%5Cmathbf%7BC%7D%5E%7B-1%7D%20%5Cmathbf%7BX%7D
"\\mathbf{X}^* = \\mathbf{C}^{-1} \\mathbf{X}") as an equivalent
estimator.

Splines, kernels, and random effects with quadratic penalties have the
form

<center>

  
![\\hat{\\beta} = (\\mathbf{X}^T\\mathbf{X} + \\lambda
\\mathbf{D}^T\\mathbf{D})^{-1}
\\mathbf{X}^T\\mathbf{Y}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cbeta%7D%20%3D%20%28%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BX%7D%20%2B%20%5Clambda%20%5Cmathbf%7BD%7D%5ET%5Cmathbf%7BD%7D%29%5E%7B-1%7D%20%5Cmathbf%7BX%7D%5ET%5Cmathbf%7BY%7D
"\\hat{\\beta} = (\\mathbf{X}^T\\mathbf{X} + \\lambda \\mathbf{D}^T\\mathbf{D})^{-1} \\mathbf{X}^T\\mathbf{Y}")  

</center>

for some smoothing parameter
![\\lambda](http://chart.apis.google.com/chart?cht=tx&chl=%5Clambda
"\\lambda"), penalty matrix
![D](http://chart.apis.google.com/chart?cht=tx&chl=D "D"). With the data
augmentation trick this can be shown to be identical to the regression
of

<center>

![\\begin{align} \\tilde{\\mathbf{Y}} = \\begin{bmatrix} y \\\\ 0
\\end{bmatrix} \\; \\; \\text{ on } \\tilde{\\mathbf{X}} =
\\begin{bmatrix} \\mathbf{X}\\\\ \\sqrt{\\lambda} \\mathbf{D}
\\end{bmatrix}
\\end{align}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbegin%7Balign%7D%20%20%5Ctilde%7B%5Cmathbf%7BY%7D%7D%20%3D%20%5Cbegin%7Bbmatrix%7D%20%20y%20%5C%5C%20%200%20%20%5Cend%7Bbmatrix%7D%20%5C%3B%20%5C%3B%20%5Ctext%7B%20on%20%7D%20%5Ctilde%7B%5Cmathbf%7BX%7D%7D%20%3D%20%5Cbegin%7Bbmatrix%7D%20%20%5Cmathbf%7BX%7D%5C%5C%20%20%5Csqrt%7B%5Clambda%7D%20%5Cmathbf%7BD%7D%20%20%5Cend%7Bbmatrix%7D%20%20%5Cend%7Balign%7D
"\\begin{align}  \\tilde{\\mathbf{Y}} = \\begin{bmatrix}  y \\\\  0  \\end{bmatrix} \\; \\; \\text{ on } \\tilde{\\mathbf{X}} = \\begin{bmatrix}  \\mathbf{X}\\\\  \\sqrt{\\lambda} \\mathbf{D}  \\end{bmatrix}  \\end{align}")

</center>

that is

<center>

![\\hat{\\beta} = (\\tilde{\\mathbf{X}}^T \\tilde{\\mathbf{X}})^T
\\tilde{\\mathbf{X}}^T
\\tilde{\\mathbf{Y}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cbeta%7D%20%3D%20%28%5Ctilde%7B%5Cmathbf%7BX%7D%7D%5ET%20%5Ctilde%7B%5Cmathbf%7BX%7D%7D%29%5ET%20%5Ctilde%7B%5Cmathbf%7BX%7D%7D%5ET%20%5Ctilde%7B%5Cmathbf%7BY%7D%7D
"\\hat{\\beta} = (\\tilde{\\mathbf{X}}^T \\tilde{\\mathbf{X}})^T \\tilde{\\mathbf{X}}^T \\tilde{\\mathbf{Y}}")

</center>

More generally, semiparametric regression models with non-normal
responses often have the form

<center>

![\\hat{\\beta} = (\\mathbf{X}^T \\mathbf{W} \\mathbf{X} + \\lambda
\\mathbf{D}^T \\mathbf{D})^T \\mathbf{X}^T \\mathbf{W}
\\mathbf{z}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cbeta%7D%20%3D%20%28%5Cmathbf%7BX%7D%5ET%20%5Cmathbf%7BW%7D%20%5Cmathbf%7BX%7D%20%2B%20%5Clambda%20%5Cmathbf%7BD%7D%5ET%20%5Cmathbf%7BD%7D%29%5ET%20%5Cmathbf%7BX%7D%5ET%20%5Cmathbf%7BW%7D%20%5Cmathbf%7Bz%7D
"\\hat{\\beta} = (\\mathbf{X}^T \\mathbf{W} \\mathbf{X} + \\lambda \\mathbf{D}^T \\mathbf{D})^T \\mathbf{X}^T \\mathbf{W} \\mathbf{z}")

</center>

6 where
![\\mathbf{z}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmathbf%7Bz%7D
"\\mathbf{z}") is the working variable. These may be put into the
desired format by first applying the data augmentation trick and then
whitening.
