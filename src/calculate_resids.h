#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector calculate_resids(NumericMatrix predictions,NumericVector response){
  NumericVector resids=response.size();
  NumericVector row_sums=calc_rowsums(predictions);
  resids=response - row_sums;
  return(resids);
}
