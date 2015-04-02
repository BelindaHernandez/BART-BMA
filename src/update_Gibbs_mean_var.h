#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List update_Gibbs_mean_var(NumericMatrix tree_table,NumericVector resids,double a,double sigma,double mu_mu,IntegerVector terminal_nodes,  List term_obs_tree){
  List update_params(2);
    
  NumericVector mu_ij(terminal_nodes.size());
  NumericVector Tj(terminal_nodes.size());
  NumericVector new_mean(terminal_nodes.size());
  NumericVector new_var(terminal_nodes.size());
  
  for(int k=0;k< terminal_nodes.size();k++){
    //update the node means    
    //find which observations are assigned to terminal node k in tree i
     
    IntegerVector term_obs=term_obs_tree[k];
    
    //get the number of observations in node k
    NumericVector temp_obs=resids[term_obs];
    mu_ij[k]=std::accumulate(temp_obs.begin(),temp_obs.end(), 0.0);
    Tj[k]=term_obs.size();
            
    new_mean[k]=(mu_ij[k]+a*mu_mu)/(Tj[k]+a);
  
    new_var[k]=(1/((1/(pow(sigma,2)))*(Tj[k]+a)));
    
    NumericVector temp;
    term_obs=temp;

  }  
  
  update_params[0]=new_mean;
  update_params[1]=new_var;
  
  return(update_params);
}
