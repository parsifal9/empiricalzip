//' Fits a zero-inflated generalized Poisson by EM 
//' 
//' @param x The data.
//' @export

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector  f12(NumericVector x, double m ){
  double tt;
  LogicalVector workspace(m);
  NumericVector workspace2(m);
  int aa = x.size();
  NumericVector f12(aa);
  
  
  // for(int i =0; i< m; i++){
  //   for(int j =0; j< m; j++){
  //     workspace(j) = x[i]==x[j]; 
  //   }
  //   f12(i)=m/sum(workspace);
  // }
  //works but very slow

  workspace2 = sort_unique(x); //even thought x is possibly already sorted -- "unique" does not return values in ascending order
  aa = x.size();
  for(int i =0; i< workspace2.size() ; i++){
    workspace = workspace2(i)==x;
    tt=  sum(workspace);
    for(int j =0; j< aa; j++){
      if( workspace(j)){
	f12(j) = m/tt;
      }
    }
  }
  
  return f12;
}


  
// f12 <- function(x, M){
//   y <- rep(0, length(x))
//   for(i in 1:length(x)){
//     y[i] <- M/length(x[which(x==x[i])])
//   }
//   return(y)
// }
