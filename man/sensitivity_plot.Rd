\name{sensitivity_plot}
\alias{sensitivity_plot}
\title{Plot how the distribution of possible effects changes with restrictions on the R-squared.}
\description{Takes the results of \link[DOPE]{DOPE} and plots how the uncertainty changes with various restrictions on the R-squared.}
\usage{
  sensitivity_plot(output,vname,adj=NULL)
}
\arguments{
  \item{output}{
    data.frame of simulations, usually from \link[DOPE]{DOPE}
  }
  \item{vname}{
    character or column index.  Which distribution would you like to see?
  }
  \item{adj}{
    numeric.  Set y-axis lower bound.
  }
}

\value{
  Returns a \link[ggplot2]{ggplot} object.
  
  The solid dots indicate how the certainty in the effect changes as you say "the world is not that deterministic," reducing the maximum allowable R-squared. 
  This reduces uncertainty relative to ignornace, indicated by the dashed line, at the cost of restrictiveness, indicated by the hollow points. 
  
  In the other direction, one might say that "the world is more deterministic than reflected in my model." 
  This is represented by the crossed points and increases uncertainty until the lower pessimistic bound on effect certainty. 
  The proportion of draws rejected by this lower thresholding is given by the filled diamonds.
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
  
  sensitivity_plot(dope,"V2")
}