\name{simfun}
\alias{simfun}
\alias{simfuncpp}
\title{Simulate an an augmented regression.}
\description{Takes a valid k by k covariance matrix, \link[DOPE]{augment}s it into a valid k+1 by k+1 covariance matrix, and returns the vector of regression coefficients and R-squared from the implicit linear regression.}
\usage{
  simfun(vcvm, buff = sqrt(.Machine$double.eps))
  simfuncpp(vcvm, buff = sqrt(.Machine$double.eps))
}
\arguments{
  \item{vcvm}{
    A valid k by k covariance matrix.  Assumes that the outcome variable is defined by the first row and column.
  }
  \item{buff}{
    numeric.  A buffer to avoid numeric positive non-definiteness.
  }
}

\value{
  A valid k + 1 length vector of regression coefficients (k-1 regressors, 1 control function, 1 R-squared).
}
\notes{
  Primarily for internal use by \link[DOPE]{DOPE}.  Currently implimented in both R and C++.  The C++ version is faster while the R version is easier for the expected user base to read and modify as needed.
}
\examples{
  set.seed(1234)
  corm <- RandomCormCPP(5)
  simfuncpp(corm)
}