//This code will take the output of BART-BMA (list of sums of trees) and will update the predicted values for the terminal node
//means and model variance

//first take one set of sum of trees:
//this will work for the training data only, this will be updated for test data as an optional parameter later on 
//(will also have to update external predict function for test predictions)
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_nodes(NumericMatrix tree_table){
  //NumericVector terminal_nodes;
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
  
  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==-1);
  term_nodes=term_nodes+1;
  return(wrap(term_nodes));
}
 
//################################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
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
//################################################################################################################################//
// [[Rcpp::export]]
NumericVector calc_rowsums(NumericMatrix predictions){
   arma::mat M1(predictions.begin(), predictions.nrow(), predictions.ncol(), false);
  arma::colvec predicted_values=sum(M1,1);
  return(wrap(predicted_values));
}
//############################################################################################//
// [[Rcpp::export]]
NumericVector calculate_resids(NumericMatrix predictions,NumericVector response){
  NumericVector resids=response.size();
  NumericVector row_sums=calc_rowsums(predictions);
  resids=response - row_sums;
  return(resids);
}

//################################################################################################################################//
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
  
  //only thing changing is mu_ij the terminal node mean
  //everything else was set outside the loop
  std::cout<<"inside mean var "<<new_mean[0]<<"\n";
  update_params[0]=new_mean;
  update_params[1]=new_var;
  
  return(update_params);
}
//################################################################################################################################//

// [[Rcpp::export]]
double update_sigma(double a1,double b,NumericVector resids,int n){
  NumericVector sq_resid=resids*resids;
  double ssr= std::accumulate(sq_resid.begin(),sq_resid.end(),0.0);
  
 // std::cout<<ssr;
  double shape=(a1+n/2);
  double rate =((ssr/2)+(1/b));
  //std::cout<<"shape"<<shape<<" "<<"rate"<<rate;
  RNGScope scope;
  double tau =R::rgamma(shape,1/rate);
  double sigma = sqrt((1/tau));
  return sigma;
//  return sigma;
}

//################################################################################################################################//
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
//################################################################################################################################//

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
      
     // std::cout<<"in overall sum trees loop "<<"\n";
      if(is<List>(s)){
      //  std::cout<<"it was a list! i is "<<i<<"\n";
        //if current set of trees contains more than one tree...usually does!
        List sum_tree=overall_sum_trees[i];
        
          List sum_tree_mat=overall_sum_mat[i];
        
        //save all info in list of list format the same as the trees.
        
        List term_nodes_trees(sum_tree.size());
        List term_obs_trees(sum_tree.size());
        NumericMatrix predictions(num_obs,sum_tree.size());
        
        for(int k =0;k<sum_tree.size();k++){
         // std::cout<<"now looping through tree "<<k<< " in sum of trees "<<i<<"\n";
          
          NumericMatrix tree_table=sum_tree[k];
        // std::cout<<"got tree table"<<"\n";
          NumericMatrix tree_mat=sum_tree_mat[k];
        //  std::cout<<"got tree mat"<<"\n";
          NumericVector term_nodes=find_term_nodes(tree_table);
        //  std::cout<<"got term nodes"<<"\n";
          term_nodes_trees[k]=term_nodes;
        //  std::cout<<"got tree info "<<"\n";
          List term_obs_tree(term_nodes.size());
          NumericVector term_preds(num_obs);
          
          for(int j=0;j<term_nodes.size();j++){
         //   std::cout<<"entered term node "<<j<<"in tree i"<<"\n";
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
       // std::cout<<"tree sum was a single tree "<<i<<"\n";
      //  List sum_tree2= overall_sum_trees[i];
      //  std::cout<<"got past here"<<"\n";
        NumericMatrix sum_tree=overall_sum_trees[i];
      //  std::cout<<"made the tree mat"<<"\n";
        NumericMatrix tree_mat=overall_sum_mat[i];
      //  std::cout<<"problem was here?"<<"\n";
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
  //std::cout<<"about to return "<<"\n";
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;
  return(ret);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix remove_curr_col(NumericMatrix predy,int i){
  arma::mat M=Rcpp::as<arma::mat>(predy);
  M.shed_col(i);
  NumericMatrix s=as<NumericMatrix>(wrap(M));
  return(s);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector get_new_mean(IntegerVector terminal_nodes,List new_mean_var){
  NumericVector node_means;
  for(int k=0;k<terminal_nodes.size();k++){
    NumericVector sd=new_mean_var[1];
    NumericVector temp_mean=new_mean_var[0];
    double new_mean= R::rnorm(temp_mean[k],sqrt(sd[k]));
    node_means.push_back(new_mean);
   
  }

  return(node_means);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    //for(int i=0;i<term_obs.size();i++){
    //  new_preds[term_obs[i]]=new_mean[k];
    //}
  }
  
  updated_preds[0]=tree_table;
  updated_preds[1]=new_preds;
  
  return(updated_preds);
}
//################################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));
  return(y_scaled);
}
//###########################################################################################################################//
// [[Rcpp::export]]
NumericVector get_original(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  return(original_y);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]

NumericVector find_internal_nodes(NumericMatrix treetable){
  NumericVector internal_nodes;

  for(int l=0;l<treetable.nrow();l++){    
    if(treetable(l,4)==1){
      internal_nodes.push_back(l+1);
    }
  }
  
  NumericVector internal_nodes_sort = clone(internal_nodes);
  std::sort(internal_nodes.begin(), internal_nodes.end());
  
  return(internal_nodes_sort);
}
//###########################################################################################################################//
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List get_tree_info_test_data(NumericMatrix test_data,NumericMatrix tree_data) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on...
  //should have an error check for this later on or at least check it's the same dimension
  
  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values
  
  //term_node_means is a vector storing the terminal node mean values
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  NumericVector internal_nodes=find_internal_nodes(tree_data);
  NumericVector terminal_nodes=find_term_nodes(tree_data);
  arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  NumericVector tree_predictions;
  
  //now for each internal node find the observations that belong to the terminal nodes
  
  NumericVector predictions(test_data.nrow());
  List term_obs(terminal_nodes.size());
  for(int i=0;i<terminal_nodes.size();i++){
    arma::mat subdata=testd;
    int curr_term=terminal_nodes[i];
    int row_index;
    int term_node=terminal_nodes[i];
    
    if(curr_term % 2==0){
      //term node is left daughter
      row_index=terminal_nodes[i];
    //  std::cout<<"term node is left daughter "<<terminal_nodes[i]<<" term mean is "<<term_node_means[i]<<"\n";
    }else{
      //term node is right daughter
          //  std::cout<<"term node is right daughter "<<terminal_nodes[i]<<" term mean is "<<term_node_means[i]<<"\n";
      row_index=terminal_nodes[i]-1;
    }
    
    //save the left and right node data into arma uvec
    
    arma::vec left_nodes=arma_tree.col(0);
    arma::vec right_nodes=arma_tree.col(1);
    arma::mat node_split_mat;    
    node_split_mat.set_size(0,3);

    while(row_index!=1){
      //for each terminal node work backwards and see if the parent node was a left or right node
      //append split info to a matrix 
      int rd=0;
      arma::uvec parent_node=arma::find(left_nodes == term_node);
      
      if(parent_node.size()==0){
        parent_node=arma::find(right_nodes == term_node);
        rd=1;
      }
      
      //want to cout parent node and append to node_split_mat
      
      node_split_mat.insert_rows(0,1);
      node_split_mat(0,0)=tree_data(parent_node[0],2);
      std::cout<<"split var is "<<node_split_mat(0,0)<<"\n";
      node_split_mat(0,1)=tree_data(parent_node[0],3);
      std::cout<<"split point is "<<node_split_mat(0,1)<<"\n";
      node_split_mat(0,2)=rd;
      
      row_index=parent_node[0]+1;
      term_node=parent_node[0]+1;
    }
    
    //once we have the split info, loop through rows and find the subset indexes for that terminal node!
    //then fill in the predicted value for that tree
    //double prediction = tree_data(term_node,5);
    arma::uvec pred_indices;
    int split= node_split_mat(0,0)-1;
    arma::vec tempvec = testd.col(split);
    double temp_split = node_split_mat(0,1);
    
    if(node_split_mat(0,2)==0){
      pred_indices = arma::find(tempvec <= temp_split);
    }else{
      pred_indices = arma::find(tempvec > temp_split);
      
    }
    
    arma::uvec temp_pred_indices;
    //arma::uvec col_indices=seq_len(testd.n_cols);
    arma::vec data_subset = testd.col(split);
    data_subset=data_subset.elem(pred_indices);
    
    //now loop through each row of node_split_mat
    int n=node_split_mat.n_rows;
    //std::cout<<"n is "<<n<<"\n";
    
    for(int j=1;j<n;j++){
      int curr_sv=node_split_mat(j,0);
      double split_p = node_split_mat(j,1);
      
      data_subset = testd.col(curr_sv-1);
      data_subset=data_subset.elem(pred_indices);
     
      if(node_split_mat(j,2)==0){
        //split is to the left
        temp_pred_indices=arma::find(data_subset <= split_p);
        std::cout<<"the split was left "<<temp_pred_indices.size()<<"\n";
      }else{
        //split is to the right
        temp_pred_indices=arma::find(data_subset > split_p);
       // std::cout<<"the spliut was right "<<temp_pred_indices.size()<<"\n";
      }
      pred_indices=pred_indices.elem(temp_pred_indices);
      
      if(pred_indices.size()==0){
        //std::cout<<"no observations for this node "<<"\n";
        continue;
      }
    
    }
   
    double nodemean=tree_data(terminal_nodes[i]-1,5);
    //std::cout<<"got the terminal node mean "<<nodemean<<"\n";
    IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
    //std::cout<<"converted predindices "<<"\n";
    //std::cout<<"predictions size is "<<predictions.size()<<"\n";
    predictions[predind]= nodemean;
   term_obs[i]=predind;
  } 
  List ret(3);
  ret[0] = terminal_nodes;
  ret[1] = term_obs;
  ret[2] = predictions;
  std::cout<<"predictions for both observations: "<<predictions[0]<<" "<<predictions[1]<<"\n";
  return(ret);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
      
     // std::cout<<"in overall sum trees loop "<<"\n";
      if(is<List>(s)){
      //  std::cout<<"it was a list! i is "<<i<<"\n";
        //if current set of trees contains more than one tree...usually does!
        List sum_tree=overall_sum_trees[i];
        
        //save all info in list of list format the same as the trees.
        
        List term_nodes_trees(sum_tree.size());
        List term_obs_trees(sum_tree.size());
        NumericMatrix predictions(num_obs,sum_tree.size());
         
        for(int k =0;k<sum_tree.size();k++){
         // std::cout<<"now looping through tree "<<k<< " in sum of trees "<<i<<"\n";
          
          NumericMatrix tree_table=sum_tree[k];
        // std::cout<<"got tree table"<<"\n";
         
        //  std::cout<<"got tree mat"<<"\n";
        //  NumericVector term_nodes=find_term_nodes(tree_table);
        //  std::cout<<"got term nodes"<<"\n";
       
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
       // std::cout<<"tree sum was a single tree "<<i<<"\n";
      //  List sum_tree2= overall_sum_trees[i];
      //  std::cout<<"got past here"<<"\n";
        NumericMatrix sum_tree=overall_sum_trees[i];
      //  std::cout<<"made the tree mat"<<"\n";
       
      //  std::cout<<"problem was here?"<<"\n";
        
        
          List tree_info=get_tree_info_test_data(test_data, sum_tree) ;
          overall_term_nodes_trees[i]=tree_info[0];
          List term_obs_tree=tree_info[1];
          NumericVector term_preds=tree_info[2];
          NumericVector predictions=term_preds;
          
        overall_term_obs_trees[i]= term_obs_tree;
        overall_predictions[i]=predictions;
      }  
  }    
  //std::cout<<"about to return "<<"\n";
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;
  return(ret);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List gibbs_sampler(List overall_sum_trees,List overall_sum_mat,NumericVector y,NumericVector BIC_weights,int num_iter,int burnin,int num_obs,int num_test_obs,double a,double sigma,double mu_mu,double nu,double lambda,List resids,NumericMatrix test_data){
  //later need to separate the sums of trees in overall list into elements of sum_trees and sum_tree_mat
  //then iterate through each individual tree.
  if(burnin>=num_iter){
    throw std::range_error("Number of iterations has to be greater than the number of burn-in samples");
  }
  std::cout<<"entered the function!"<<overall_sum_trees.size() <<"\n";
  //get tree info just once don't repeat for each iteration:
  //need terminal nodes, terminal node means and observations terminal nodes refer to
  List tree_info=get_tree_info(overall_sum_trees,overall_sum_mat,num_obs);
  List test_tree_info=get_tree_info_testdata_overall(overall_sum_trees,num_test_obs,test_data);
  List overall_term_test_obs_trees=test_tree_info[1];

  double a1=nu/2;
  double b=2/(nu*lambda);
  std::cout<<"got tree info"<<"\n"; 
  double sigma_init=sigma;
  List overall_term_nodes_trees=tree_info[0];
  List overall_term_obs_trees=tree_info[1];
  List overall_predictions=tree_info[2];
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  std::cout<<"got tree info2"<<"\n"; 
  List prediction_list(overall_sum_trees.size());
  List prediction_list_orig(overall_sum_trees.size());
  List prediction_test_list(overall_sum_trees.size());
  List prediction_test_list_orig(overall_sum_trees.size());
  List overall_sum_trees1=clone(overall_sum_trees);
  List overall_sum_mat1=clone(overall_sum_mat);
  std::cout<<"got tree info3 "<<overall_sum_trees.size()<<"\n"; 
  List sigma_chains(overall_sum_trees.size());
   std::cout<<"got tree info4"<<"\n"; 
//  IntegerVector tn=overall_term_nodes_trees[0];
   std::cout<<"got tree info5"<<"\n"; 
//  NumericMatrix muijtree1(num_iter,tn.size()); 
  std::cout<<"about to loop through sum of trees"<<"\n";
  int one_tree=0;
 // List muijs(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //for each set of trees loop over individual trees    
    NumericVector sigma_its(num_iter);
    sigma=sigma_init;
    SEXP s = overall_sum_trees[i];
    NumericVector test_preds_sum_tree;
    NumericMatrix sum_predictions;
    NumericMatrix sum_test_predictions;
    std::cout<<"got sum of trees for list "<<i<<"\n";
    
    
      
      if(is<List>(s)){
        std::cout<<"now going to update node means"<<"\n";
        //if current set of trees contains more than one tree...usually does!
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
       std::cout<<"matricews initialised"<<"\n";
        NumericMatrix sum_new_predictions(sum_predictions.nrow(),sum_predictions.ncol());
         NumericMatrix sum_new_test_predictions(sum_predictions.nrow(),sum_predictions.ncol());
    //    post_predictions=temp;
    //    post_predictions_orig=temp;
        for(int j=0;j<num_iter;j++){
          for(int k =0;k<sum_tree.size();k++){
          std::cout<<"tree "<<j<<"is list of trees looping through "<<k<<"th tree"<<"\n";
          NumericMatrix tree_table=sum_tree[k];
          IntegerMatrix tree_mat=sum_tree_mat[k];
          //find terminal node means and observations associated with them
          IntegerVector term_nodes=sum_term_nodes[k];
          List term_obs=sum_term_obs[k];
           List term_test_obs=sum_term_test_obs[k];
        //  NumericVector predictions=sum_predictions(_,k);
          NumericVector predictions=sum_resids[k];
          //current predictions are the residuals for sum of trees!
            
          //update the means and predictions for tree
          List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma,mu_mu,term_nodes,term_obs);
        //  NumericVector new_node_mean=new_node_mean_var[0];
           NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);
          NumericVector new_node_var=new_node_mean_var[1];
          std::cout<<"updated mean and var"<<"\n";
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          
          List updated_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
          NumericVector temp_preds=updated_preds[1];
          sum_new_predictions(_,k)=temp_preds;
        
         
          List updated_test_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_test_obs);         
          NumericVector temp_test_preds=updated_test_preds[1];
          sum_new_test_predictions(_,k)=temp_test_preds;
          std::cout<<"and predictions..."<<"\n";
          //update sigma
//          NumericVector pred_obs1=calc_rowsums(sum_predictions);
         // NumericVector S=calculate_resids(sum_predictions,y_scaled);  
          
         /* if(i==0&&j==0){
            List ret(3);
            ret[0]=sum_predictions;
            ret[1]=S;
            ret[2]=y_scaled;
            return(ret);
          }
        */  
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

        std::cout<<"left iterations of sum of trees "<<i<<" post pred dims "<<post_predictions.nrow()<<" "<<post_predictions.ncol() <<"\n";
        prediction_list[i]=post_predictions;
        prediction_list_orig[i]=post_predictions_orig;
        prediction_test_list[i]=post_test_predictions;
        prediction_test_list_orig[i]=post_test_predictions_orig;
        std::cout<<"set predictions"<<"\n";
      }
     sigma_chains[i]=sigma_its;
      
      }else{
        
        one_tree=1;
      }  
      
      
    
    
   /* if(i==0){
      List ret(2);
      ret[0]= prediction_list;
      ret[1]= prediction_list_orig;
      return(ret);
    }
  */  
    
  }
  
  if(one_tree==1){
    
    std::cout<<"only one tree!"<<"\n";
    NumericVector sigma_its(num_iter);
    NumericMatrix post_predictions(num_iter,num_obs);
    NumericMatrix post_predictions_orig(num_iter,num_obs);
    NumericMatrix post_test_predictions(num_iter,num_obs);
    NumericMatrix post_test_predictions_orig(num_iter,num_obs);

    std::cout<<"overall predictions size is "<<overall_predictions.size()<<" tree size "<<overall_sum_trees.size()<<"\n";
    NumericMatrix sum_predictions(num_obs,overall_predictions.size());
     NumericMatrix sum_test_predictions(num_test_obs,overall_predictions.size());
    for(int t=0;t<overall_predictions.size();t++){
          NumericVector preds=overall_predictions[t];
          sum_predictions(_,t)=preds;
    }
 //   NumericVector overallsigma(num_iter*overall_sum_trees.size());
 //   int count=0;
    for(int j=0;j<num_iter;j++){
    for(int i=0;i<overall_sum_trees.size();i++){
      
      NumericMatrix tree_table=overall_sum_trees1[i];
      IntegerMatrix tree_mat=overall_sum_mat1[i];
      IntegerVector term_nodes=overall_term_nodes_trees[i];
      List term_obs=overall_term_obs_trees[i];
      List term_test_obs=overall_term_test_obs_trees[i];
     // NumericVector predictions=sum_predictions(_,i);
     std::cout<<"about to get predictions";
      NumericVector predictions=resids[i];
      std::cout<<"got them";
   //   NumericMatrix temp_predy=remove_curr_col(sum_predictions,i);
      
  //    NumericVector predictions=calculate_resids(temp_predy,y_scaled);
      
      
          //find terminal node means and observations associated with them
          
          //current predictions are the residuals for sum of trees!
            
          //update the means and predictions for tree
          List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma,mu_mu,term_nodes,term_obs);
           NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);
        //  NumericVector new_node_mean=new_node_mean_var[0];
          NumericVector new_node_var=new_node_mean_var[1];
          //if(i==0){
         //   muijtree1(j,_)=new_node_mean;
         // }
          std::cout<<"updated mean and var "<<new_node_mean[0]<<" "<<sigma<<"\n";
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          std::cout<<"pred before update "<< sum_predictions(0,i)<<"\n";
          List updated_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
          NumericVector temp_preds=updated_preds[1];
          sum_predictions(_,i)=temp_preds;
        
        //get updated predictions for the test data
         List updated_test_preds=update_predictions(tree_table,new_node_mean,new_node_var,num_test_obs,term_nodes,term_test_obs);         
        NumericVector temp_test_preds=updated_test_preds[1];
        sum_predictions(_,i)=temp_preds;
         sum_test_predictions(_,i)=temp_test_preds;
        //  std::cout<<"and predictions..."<<temp_preds[0]<<"\n";
          //update sigma
        //NumericVector pred_obs1=calc_rowsums(sum_predictions);
          NumericVector S=calculate_resids(sum_predictions,y_scaled);  
      
          //get overall predictions for current iteration and current sum of trees
          sigma= update_sigma(a1,b,S,num_obs);
          sigma_its[j]=sigma;
      //  overallsigma[count]=sigma;
     //   count++;
        
       
    //    std::cout<<"left iterations of sum of trees "<<i<<" post pred dims "<<post_predictions.nrow()<<" "<<post_predictions.ncol() <<"\n";
        
  //      std::cout<<"set predictions"<<"\n";
       // muijs[i]=muijtree;
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
 // sigma_chains[0]=overallsigma;
  std::cout<<"sigma its "<<sigma_its[0]<<" "<<sigma_its[1]<<"\n";
  NumericVector test2=sigma_chains[0];
  std::cout<<"sigma chains "<<test2[0]<<" "<<test2[1]<<"\n";
  prediction_list[0]=post_predictions;
  std::cout<<"got here"<<"\n";
  prediction_list_orig[0]=post_predictions_orig;
   prediction_test_list[0]=post_test_predictions;
  prediction_test_list_orig[0]=post_test_predictions_orig;
   std::cout<<"about to go around again!"<<"\n";
  }
  
  List ret(5);
  NumericVector test2=sigma_chains[0];
  std::cout<<"sigma chains1 "<<test2[0]<<" "<<test2[1]<<"\n";
  ret[0]= prediction_list;
  ret[1]= prediction_list_orig;
  ret[2]=sigma_chains;
  ret[3]=prediction_test_list;
  ret[4]=prediction_test_list_orig;
//  ret[3]=muijtree1;
  return(ret);
  
  std::cout<<"at very end going to return"<<"\n";
  //at the end add all the predicted values up to get the overall prediction 
 
}
