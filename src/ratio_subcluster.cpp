#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cassert>
#include <functional>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double subcluster_stable(const arma::uvec& ind, const arma::mat& M, const arma::mat& N, const int steps) {

  // ARGUMENTS

  int c = ind.size();
  int n = M.n_rows;
  arma::rowvec V(n, arma::fill::zeros);

  V(ind-1) = arma::ones<arma::rowvec>(c);

  arma::rowvec out(n);
  arma::rowvec out_temp(n);
  out = V;

  for (int j = 0; j < steps; j++) {
    out_temp = out;
    for (int i = 0; i < n; i++) {
      arma::rowvec x = M.row(i);
      arma::rowvec y = N.row(i);
      double xv = arma::sum(out_temp.elem(arma::conv_to<arma::uvec>::from(x - 1)) % y.t());
      out(i) = xv/arma::sum(y);
    }
  }

  Rcpp::NumericVector out_vec = Rcpp::NumericVector(out.begin(), out.end());
  std::nth_element(out_vec.begin(), out_vec.begin()+c-1, out_vec.end(), std::greater<double>());
  double cut = out_vec(c-1);
  arma::uvec idd = ind(arma::find(out(ind-1) >= cut));
  double a = (double)idd.size();
  double b = (double)ind.size();
  double xx = a/b;
  // return (Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(out)));
  return (xx);
}
