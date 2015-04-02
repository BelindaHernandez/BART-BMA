#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_tree_info_testdata_overall(List overall_sum_trees,int num_obs,NumericMatrix test_data){
    
  List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  List overall_predictions(overall_sum_trees.size());
        
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      List sum_tree=overall_sum_trees[i];        
      //save all info in list of list format the same as the trees.       
      List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      NumericMatrix predictions(num_obs,sum_tree.size());
       
      for(int k =0;k<sum_tree.size();k++){          
        NumericMatrix tree_table=sum_tree[k];
        List tree_info=get_tree_info_test_data(test_data, tree_table) ;
        NumericVector term_nodes=tree_info[0];
        term_nodes_trees[k]=term_nodes;
        term_obs_trees[k]=tree_info[1];
        NumericVector term_preds=tree_info[2];
        predictions(_,k)=term_preds;
      } 
      overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      List tree_info=get_tree_info_test_data(test_data, sum_tree) ;
      overall_term_nodes_trees[i]=tree_info[0];
      List term_obs_tree=tree_info[1];
      NumericVector term_preds=tree_info[2];
      NumericVector predictions=term_preds;   
      overall_term_obs_trees[i]= term_obs_tree;
      overall_predictions[i]=predictions;
    }  
  }    
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;
  return(ret);
}
