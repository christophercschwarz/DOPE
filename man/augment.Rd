\name{augment}
\alias{augment}
\alias{augmentcpp}
\title{Augment a covariance matrix.}
\description{Takes a valid $k \times k$ covariance matrix and produces a valid $k+1 \times k+1$ covariance matrix, uniformly drawn from that set, with the upper $k \times k$ submatrix being the input.}
\usage{
augment(covm, buff = sqrt(.Machine$double.eps))
augmentcpp(covm, buff = sqrt(.Machine$double.eps))
}
\arguments{
  \item{covm}{
    A valid $k \times k$ covariance matrix.
  }
  \item{buff}{
    numeric.  A buffer to avoid numeric positive non-definiteness.
  }
}

\value{
   A valid $k \times k$ covariance matrix.
}
\notes{
  Currently implimented in both R and C++.  The C++ version is faster while the R version is easier for the expected user base to read and modify as needed.
}
\examples{
  set.seed(1234)
  corm <- RandomCormCPP(5)
  augmentcpp(corm)
}

