#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec lrm_coef(const mat& X, const vec& Y){
  vec coef=solve(X, Y);
  return(coef);
}





// test case
// X.vec <- mtcars$wt
// X.mat <- cbind(rep(1,length(X.vec)), X.vec)
// Y <- mtcars$mpg
// lrm_coef(X.mat, Y)

