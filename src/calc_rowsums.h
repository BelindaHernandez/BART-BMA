#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector calc_rowsums(NumericMatrix predictions){
  arma::mat M1(predictions.begin(), predictions.nrow(), predictions.ncol(), false);
  arma::colvec predicted_values=sum(M1,1);
  return(wrap(predicted_values));
}
