#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List gibbs_sampler(List overall_sum_trees,List overall_sum_mat,NumericVector y,NumericVector BIC_weights,int num_iter,int burnin,int num_obs,int num_test_obs,double a,double sigma,double mu_mu,double nu,double lambda,List resids,NumericMatrix test_data){
  if(burnin>=num_iter){
    throw std::range_error("Number of iterations has to be greater than the number of burn-in samples");
  }
  //get tree info just once don't repeat for each iteration:
  //need terminal nodes, terminal node means and observations terminal nodes refer to
  List tree_info=get_tree_info(overall_sum_trees,overall_sum_mat,num_obs);
  List test_tree_info=get_tree_info_testdata_overall(overall_sum_trees,num_test_obs,test_data);
  List overall_term_test_obs_trees=test_tree_info[1];

  double a1=nu/2;
  double b=2/(nu*lambda);
  double sigma_init=sigma;
  List overall_term_nodes_trees=tree_info[0];
  List overall_term_obs_trees=tree_info[1];
  List overall_predictions=tree_info[2];
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  List prediction_list(overall_sum_trees.size());
  List prediction_list_orig(overall_sum_trees.size());
  List prediction_test_list(overall_sum_trees.size());
  List prediction_test_list_orig(overall_sum_trees.size());
  List overall_sum_trees1=clone(overall_sum_trees);
  List overall_sum_mat1=clone(overall_sum_mat); 
  List sigma_chains(overall_sum_trees.size());
  int one_tree=0;
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //for each set of trees loop over individual trees    
    NumericVector sigma_its(num_iter);
    sigma=sigma_init;
    SEXP s = overall_sum_trees[i];
    NumericVector test_preds_sum_tree;
    NumericMatrix sum_predictions;
    NumericMatrix sum_test_predictions;
   
      if(is<List>(s)){
        //if current set of trees contains more than one tree
        List sum_tree=overall_sum_trees1[i];
        List sum_tree_mat=overall_sum_mat1[i];
        List sum_term_nodes=overall_term_nodes_trees[i];
        List sum_term_obs=overall_term_obs_trees[i];
        List sum_term_test_obs=overall_term_test_obs_trees[i];
        List sum_resids=resids[i];
        NumericMatrix tree_predictions=overall_predictions[i];
        sum_predictions=clone(tree_predictions);
        NumericMatrix post_predictions(num_iter,num_obs);
        NumericMatrix post_predictions_orig(num_iter,num_obs);
        NumericMatrix post_test_predictions(num_iter,num_test_obs);
        NumericMatrix post_test_predictions_orig(num_iter,num_test_obs);
        NumericMatrix sum_new_predictions(sum_predictions.nrow(),sum_predictions.ncol());
        NumericMatrix sum_new_test_predictions(sum_predictions.nrow(),sum_predictions.ncol());

        for(int j=0;j<num_iter;j++){
          for(int k =0;k<sum_tree.size();k++){
          NumericMatrix tree_table=sum_tree[k];
          IntegerMatrix tree_mat=sum_tree_mat[k];
          //find terminal node means and observations associated with them
          IntegerVector term_nodes=sum_term_nodes[k];
          List term_obs=sum_term_obs[k];
          List term_test_obs=sum_term_test_obs[k];
          NumericVector predictions=sum_resids[k];
          //current predictions are the residuals for sum of trees!
            
          //update the means and predictions for tree
          List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma,mu_mu,term_nodes,term_obs);
          NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);
          NumericVector new_node_var=new_node_mean_var[1];
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          
          List updated_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
          NumericVector temp_preds=updated_preds[1];
          sum_new_predictions(_,k)=temp_preds;
        
         
          List updated_test_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_test_obs);         
          NumericVector temp_test_preds=updated_test_preds[1];
          sum_new_test_predictions(_,k)=temp_test_preds;
 
          //get overall predictions for current iteration and current sum of trees
          sigma= update_sigma(a1,b,predictions,num_obs);
          sigma_its[j]=sigma;
        }
        NumericVector pred_obs=calc_rowsums(sum_new_predictions);
        post_predictions(j,_)=pred_obs;
        NumericVector original_y=get_original(min(y),max(y),-0.5,0.5,pred_obs);    
        post_predictions_orig(j,_)=original_y;

        NumericVector pred_test_obs=calc_rowsums(sum_new_test_predictions);
        post_test_predictions(j,_)=pred_test_obs;
        NumericVector original_test_y=get_original(min(y),max(y),-0.5,0.5,pred_test_obs);    
        post_test_predictions_orig(j,_)=original_test_y;

        prediction_list[i]=post_predictions;
        prediction_list_orig[i]=post_predictions_orig;
        prediction_test_list[i]=post_test_predictions;
        prediction_test_list_orig[i]=post_test_predictions_orig;
      }
      sigma_chains[i]=sigma_its;
      
      }else{        
        one_tree=1;
      }      
  }
  
  if(one_tree==1){
    NumericVector sigma_its(num_iter);
    NumericMatrix post_predictions(num_iter,num_obs);
    NumericMatrix post_predictions_orig(num_iter,num_obs);
    NumericMatrix post_test_predictions(num_iter,num_obs);
    NumericMatrix post_test_predictions_orig(num_iter,num_obs);
    NumericMatrix sum_predictions(num_obs,overall_predictions.size());
    NumericMatrix sum_test_predictions(num_test_obs,overall_predictions.size());
    for(int t=0;t<overall_predictions.size();t++){
      NumericVector preds=overall_predictions[t];
      sum_predictions(_,t)=preds;
    }

    for(int j=0;j<num_iter;j++){
      for(int i=0;i<overall_sum_trees.size();i++){     
        NumericMatrix tree_table=overall_sum_trees1[i];
        IntegerMatrix tree_mat=overall_sum_mat1[i];
        IntegerVector term_nodes=overall_term_nodes_trees[i];
        List term_obs=overall_term_obs_trees[i];
        List term_test_obs=overall_term_test_obs_trees[i];
        NumericVector predictions=resids[i];
        //find terminal node means and observations associated with them
            
        //current predictions are the residuals for sum of trees
              
        //update the means and predictions for tree
        List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma,mu_mu,term_nodes,term_obs);
        NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);      
        NumericVector new_node_var=new_node_mean_var[1];
        List updated_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
        NumericVector temp_preds=updated_preds[1];
        sum_predictions(_,i)=temp_preds;
          
        //get updated predictions for the test data
        List updated_test_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_test_obs,term_nodes,term_test_obs);         
        NumericVector temp_test_preds=updated_test_preds[1];
        sum_predictions(_,i)=temp_preds;
        sum_test_predictions(_,i)=temp_test_preds;
        NumericVector S=calculate_resids(sum_predictions,y_scaled);  
        //get overall predictions for current iteration and current sum of trees
        sigma= update_sigma(a1,b,S,num_obs);
        sigma_its[j]=sigma;
      }
      NumericVector pred_obs=calc_rowsums(sum_predictions);
      NumericVector pred_test_obs=calc_rowsums(sum_test_predictions);
      post_predictions(j,_)=pred_obs;
      post_test_predictions(j,_)=pred_test_obs;
      NumericVector original_y=get_original(min(y),max(y),-0.5,0.5,pred_obs);   
      NumericVector original_test_y=get_original(min(y),max(y),-0.5,0.5,pred_test_obs);   
      post_predictions_orig(j,_)=original_y;
      post_test_predictions_orig(j,_)=original_test_y;
    }
  
    List g(1);
    sigma_chains=clone(g);
    prediction_list=clone(g);
    prediction_list_orig=clone(g);
    prediction_test_list=clone(g);
    prediction_test_list_orig=clone(g);
    sigma_chains[0]=sigma_its;
    NumericVector test2=sigma_chains[0];
    prediction_list[0]=post_predictions;
    prediction_list_orig[0]=post_predictions_orig;
    prediction_test_list[0]=post_test_predictions;
    prediction_test_list_orig[0]=post_test_predictions_orig;
  }
  
  List ret(5);
  NumericVector test2=sigma_chains[0];
  ret[0]= prediction_list;
  ret[1]= prediction_list_orig;
  ret[2]=sigma_chains;
  ret[3]=prediction_test_list;
  ret[4]=prediction_test_list_orig;
  return(ret); 
}
