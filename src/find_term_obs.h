#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector find_term_obs(NumericMatrix tree_matrix_temp,double terminal_node){
  //NumericVector term_obs;
  arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
  arma::uvec term_obs;
  
  for(int j=0;j<tree_matrix_temp.ncol();j++){
    arma::vec colmat=arma_tree_mat.col(j);
    term_obs=arma::find(colmat==terminal_node);
    if(term_obs.size()>0){
      break;
    }
  }
  
  return(wrap(term_obs));
}
