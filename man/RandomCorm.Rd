\name{RandomCorm}
\alias{RandomCorm}
\alias{RandomCormCPP}
\title{Generate a valid correlation matrix}
\description{Generates a valid correlation matrix of desired size.}
\usage{
  RandomCorm(nvars, buff = sqrt(.Machine$double.eps))
  RandomCormCPP(nvars, buff = sqrt(.Machine$double.eps))
}
\arguments{
  \item{nvars}{
    numeric size of desired correlation matrix
  }
  \item{buff}{
    numeric.  A buffer to avoid numeric positive non-definiteness.
  }
}

\value{
  A valid $nvars \times nvars$ correlation matrix.
}
\notes{
  Repeatedly \link[DOPE]{augment}s until desired size is met.  Currently implimented in both R and C++.  The C++ version is faster while the R version is easier for the expected user base to read and modify as needed.
}
\examples{
  set.seed(1234)
  RandomCormCPP(5)
}