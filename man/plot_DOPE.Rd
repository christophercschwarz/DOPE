\name{plot_DOPE}
\alias{plot_DOPE}
\title{Plot a distribution of possible effects.}
\description{Takes the results of \link[DOPE]{DOPE} and plots the distributions with usually good defaults.}
\usage{
  plot_DOPE(output,vname,xmin=NULL,xmax=NULL,bw=NULL,shade=FALSE,include_naive = TRUE)
}
\arguments{
  \item{output}{
    data.frame of simulations, usually from \link[DOPE]{DOPE}
  }
  \item{vname}{
    character or column index.  Which distribution would you like to see?
  }
  \item{xmin}{
    numeric.  Change the lower bound of plotted distribution.
  }
  \item{xmax}{
    numeric.  Change the upper bound of plotted distribution.
  }
  \item{bw}{
    numeric.  Set distribution bin width.
  }
  \item{shade}{
    logical.  Shade the distribution according to R-squared values?
  }
  \item{include_naive}{
    logical.  If it exists, indicate the naive regression coefficient?
  }
}

\value{
  Returns a \link[ggplot2]{ggplot} object of the simulated distribution of possible effects with numerical summaries.
}
\examples{
  set.seed(8675309)
  x_vars <- 5
  n_obs <- 1000
  corm <- RandomCormCPP(nvars = x_vars)
  X_mat <- MASS::mvrnorm(n_obs, rep(0,x_vars), Sigma = corm, empirical = TRUE)
  
  betas <- 1:x_vars
  
  y <- X_mat %*% betas + rnorm(n_obs, 0, 1)
  
  dat <- data.frame(y,X_mat)
  
  mod <- lm(y ~ ., data=dat)
  dope <- DOPE(mod, nsims = 3000, n.cores = parallel::detectCores())
  
  plot_DOPE(dope,"V2")
}