#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
// [[Rcpp::depends(RcppArmadillo)]]

NumericVector find_term_nodes(NumericMatrix tree_table){
  //NumericVector terminal_nodes;
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
  
  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==-1);
  term_nodes=term_nodes+1;
  return(wrap(term_nodes));
}
