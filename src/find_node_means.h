#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_node_means(NumericMatrix sum_tree, NumericVector term_nodes){
  NumericVector means=sum_tree(_,5);
  arma::vec node_means = as<arma::vec>(means);
  arma::uvec arma_term_nodes=as<arma::uvec>(term_nodes);
  arma_term_nodes=arma_term_nodes-1;
  arma::vec final_means=node_means.elem(arma_term_nodes);
  return(wrap(final_means));
  
}
