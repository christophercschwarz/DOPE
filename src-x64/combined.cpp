#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export()]]
arma::mat vecmult (arma::mat x, arma::mat y) {
  for(unsigned int j=0; j < x.size(); j++){
    x[j] = x[j] * y[j];
  }
  return(x);
}

// [[Rcpp::export()]]
arma::mat subset (arma::mat x, int row, int column_min, int column_max) {
  arma::mat submat = x.row(row-1);
  NumericVector out(column_max-column_min+1);
  for(int j = (column_min-1); j < column_max; j++){
    out[j] = submat[j];
  }
  return(out);
}

// [[Rcpp::export()]]
arma::mat matprod (arma::mat x) {
  arma::mat out = x * x.t();
  return(out);
}

// [[Rcpp::export()]]
arma::mat augment_loop (arma::mat x, int n, double buff) {
  arma::mat tmp = x;
  int rn = n+1;
  for (int j=1; j<n; j++){
    arma::mat B_up = x;
    B_up(n,n) = 0;
    arma::mat B_lo = x;
    B_lo(n,n) = 0;
    arma::mat cholrow = subset(x,rn,1,j);
    double sumsq = std::inner_product(cholrow.begin(),cholrow.end(),cholrow.begin(),0.0);
    B_up(n,j) = sqrt(1-sumsq);
    B_lo(n,j) = -sqrt(1-sumsq);
    arma::mat C_up = B_up * B_up.t();
    double upper = C_up(n,j);
    arma::mat C_lo = B_lo * B_lo.t();
    double lower = C_lo(n,j);
    double diff = upper - lower - buff;
    double cor = 80;
    if (diff <= 0)
        cor = (upper + lower)/2;
    if (diff > 0)
        cor = R::runif(lower,upper);
    arma::vec vect = subset(x,j+1,1,j);
    arma::mat tmp = vecmult(cholrow,vect);
    x(n,j) = 1/x(j,j) * (cor - accu(tmp));
    arma::mat cholrow2 = subset(x,rn,1,n);
    double sumsq2 = std::inner_product(cholrow2.begin(),cholrow2.end(),cholrow2.begin(),0.0);
    x(n,n) = sqrt(1 - sumsq2);
  }
  return x;
}
