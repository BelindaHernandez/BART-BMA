#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_tree_info(List overall_sum_trees,List overall_sum_mat,int num_obs){
  List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  List overall_predictions(overall_sum_trees.size());
        
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
  
    NumericVector test_preds_sum_tree;
    
    if(is<List>(s)){
      //if current set of trees contains more than one tree
      List sum_tree=overall_sum_trees[i];
      
      List sum_tree_mat=overall_sum_mat[i];
      
      //save all info in list of list format the same as the trees.
      
      List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k =0;k<sum_tree.size();k++){
        
        NumericMatrix tree_table=sum_tree[k];
        NumericMatrix tree_mat=sum_tree_mat[k];
        NumericVector term_nodes=find_term_nodes(tree_table);
        term_nodes_trees[k]=term_nodes;
        List term_obs_tree(term_nodes.size());
        NumericVector term_preds(num_obs);
        
        for(int j=0;j<term_nodes.size();j++){
          double terminal_node= term_nodes[j]; 
          NumericVector term_obs=find_term_obs(tree_mat,terminal_node);
          NumericVector node_means=find_node_means(tree_table,term_nodes);
          term_obs_tree[j]=term_obs;
          double node_mean=node_means[j];
          term_preds[term_obs]=node_mean; 
        }          
        term_obs_trees[k]=term_obs_tree;
                   
        predictions(_,k)=term_preds;
      } 
      overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      NumericMatrix tree_mat=overall_sum_mat[i];
      NumericVector term_nodes=find_term_nodes(sum_tree);
      NumericVector node_means=find_node_means(sum_tree,term_nodes);
      List term_obs_tree(term_nodes.size());
      overall_term_nodes_trees[i]=term_nodes;
      NumericVector predictions(num_obs);
      
      for(int j=0;j<term_nodes.size();j++){
        double terminal_node= term_nodes[j];
        double node_mean=node_means[j];
        NumericVector term_obs=find_term_obs(tree_mat,terminal_node);
        term_obs_tree[j]=term_obs;
        predictions[term_obs]=node_mean;
      }
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
