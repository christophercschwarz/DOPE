\name{DOPE}
\alias{DOPE}
\title{Generate a distribution of possible effects.}
\description{Wrapper for parallel \link[DOPE]{simfun}.  Takes a linear regression model fit by \link[stats]{lm} and returns the results information on the distribution of possible effects.  Currently implimented in both R and C++.  The C++ version is faster while the R version is easier for the expected user base to read and modify as needed.}
\usage{
  DOPE(mod,nsims=10000,language="cpp",n.cores=1,buff=sqrt(.Machine$double.eps))
}
\arguments{
  \item{mod}{
    A linear regression model fit by \link[stats]{lm}.
  }
  \item{nsims}{
    numeric.  How many draws to take?
  }
  \item{language}{
    character, either "cpp" or "R" determining which implimentation to use.
  }
  \item{n.cores}{
    numeric.  How many cores should the simulation be run on?
  }
  \item{buff}{
    numeric.  A buffer to avoid numeric positive non-definiteness.
  }
}

\value{
  Returns a data.frame of nsims + 1 rows, with the last row containing the input coefficients.  Includes intercept, regressor coefficients, control function coefficient, and R-squared as columns.
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
}