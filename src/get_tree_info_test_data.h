#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List get_tree_info_test_data(NumericMatrix test_data,NumericMatrix tree_data) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on
  
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
    }else{
      //term node is right daughter
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
      node_split_mat(0,1)=tree_data(parent_node[0],3);
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
    arma::vec data_subset = testd.col(split);
    data_subset=data_subset.elem(pred_indices);
    
    //now loop through each row of node_split_mat
    int n=node_split_mat.n_rows;
    
    for(int j=1;j<n;j++){
      int curr_sv=node_split_mat(j,0);
      double split_p = node_split_mat(j,1);
      
      data_subset = testd.col(curr_sv-1);
      data_subset=data_subset.elem(pred_indices);
     
      if(node_split_mat(j,2)==0){
        //split is to the left
        temp_pred_indices=arma::find(data_subset <= split_p);
      }else{
        //split is to the right
        temp_pred_indices=arma::find(data_subset > split_p);
      }
      pred_indices=pred_indices.elem(temp_pred_indices);
      
      if(pred_indices.size()==0){
        continue;
      }
    
    }
   
  double nodemean=tree_data(terminal_nodes[i]-1,5);
  IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
  predictions[predind]= nodemean;
  term_obs[i]=predind;
  } 
  List ret(3);
  ret[0] = terminal_nodes;
  ret[1] = term_obs;
  ret[2] = predictions;
  return(ret);
}
