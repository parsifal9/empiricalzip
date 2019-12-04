//' Fits a zero-inflated generalized Poisson by EM 
//' 
//' @param x The data.
//' @export

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector  gpdMixture(NumericVector x,double eta, double lambda, double theta){
  double m = x.size();
  NumericVector gpd(m);
  NumericVector r(m);
  
  gpd = exp(log(lambda) + (x - 1)*log(lambda + theta*x) - lfactorial(x)  - lambda - theta*x) + pow(10,-10);
  
  r = 1.0*(x==0) ;
  r = eta*r  +(1-eta)*gpd + pow(10,-10);
  return r;
}

// gpd <- function(x, lambda, theta){
//   r <- exp(log(lambda) + (x - 1)*log(lambda + theta*x) - lfactorial(x) - lambda - theta*x) + 10^(-10)
//   return(r)
// }

// gpdMixture <- function(x,eta,lambda, theta){
//   r <- eta*(1*(x==0))+(1-eta)*gpd(x,lambda,theta) + 10^(-10)
//   return(r) 
// }
//  gpdMixture calls gpd but we have incorporated it into the function
