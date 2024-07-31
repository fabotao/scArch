#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::NumericVector MeanShift1(const arma::rowvec& V, const arma::mat& M, const arma::mat& N, const int steps) {

  // ARGUMENTS

  int c = V.size();
  arma::rowvec out(c);
  arma::rowvec out_temp(c);
  out = V;

  for (int j = 0; j < steps; j++) {
    out_temp = out;
    for (int i = 0; i < c; i++) {
      arma::rowvec x = M.row(i);

      arma::rowvec y = N.row(i);
      double xv = arma::sum(out_temp.elem(arma::conv_to<arma::uvec>::from(x - 1)) % y.t());
      out(i) = xv/arma::sum(y);
    }
    if(out.max() - out.min() < 0.00001) {
      break;
    }

  }

  return (Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(out)));

}





