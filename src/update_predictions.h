#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List update_predictions(NumericMatrix tree_table,NumericVector new_mean,NumericVector new_var,int n,
  IntegerVector terminal_nodes,List term_obs_tree){
    
  List updated_preds(2);
  NumericVector new_preds(n);
  
  for(int k=0;k<terminal_nodes.size();k++){
    
    NumericVector term_obs=term_obs_tree[k];        
    //update the mean of the selected tree nodes:
    tree_table(terminal_nodes[k]-1,5)= new_mean[k];
    tree_table(terminal_nodes[k]-1,6)=sqrt(new_var[k]);
    double newmean=new_mean[k];
    //update residuals for next iteration
    new_preds[term_obs]=newmean;
  }
  
  updated_preds[0]=tree_table;
  updated_preds[1]=new_preds;
  
  return(updated_preds);
}
