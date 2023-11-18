#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec lrm_coef(const mat& X, const vec& Y){
  vec coef=solve(X, Y);
  return(coef);
}
