#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double update_sigma(double a1,double b,NumericVector resids,int n){
  NumericVector sq_resid=resids*resids;
  double ssr= std::accumulate(sq_resid.begin(),sq_resid.end(),0.0);
  double shape=(a1+n/2);
  double rate =((ssr/2)+(1/b));
  RNGScope scope;
  double tau =R::rgamma(shape,1/rate);
  double sigma = sqrt((1/tau));
  return sigma;
}
