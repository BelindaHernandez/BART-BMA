#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
// [[Rcpp::export]]
IntegerVector csample_num( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
) {
	RNGScope scope;
	IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
	return ret;
}

//######################################################################################################################//
  
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix add_rows(NumericMatrix prior_tree_table_temp,int grow_node){
	arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
	M(grow_node-1,5)=0;
	M(grow_node-1,6)=0;
	M(grow_node-1,0)=grow_node+1;
	M(grow_node-1,1)=grow_node+2;
	M.insert_rows(grow_node,2);
	M(grow_node,4)=-1;
	M(grow_node+1,4)=-1;
	NumericMatrix t=as<NumericMatrix>(wrap(M));
	IntegerVector rname=seq_len(t.nrow());

	List dimnms = // two vec. with static names
	List::create(rname,
			   CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
	// and assign it
	t.attr("dimnames") = dimnms;

	return(t);
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix addcol(NumericMatrix prior_tree_matrix_temp,int grow_node,NumericVector ld_obs,NumericVector rd_obs){
	int ncol=prior_tree_matrix_temp.ncol();
	arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
	M.insert_cols(ncol,1);
	for(int i =0;i<ld_obs.size();i++){
		try{
		  if(ld_obs[i]>prior_tree_matrix_temp.nrow()){
			throw std::range_error("can't add col because ld row index is out of range");
		  }
		}catch(...){
		  ::Rf_error("there is a problem adding col to mat don't know why");
		}
		M(ld_obs[i],ncol)=grow_node+1;
	}
	for(int i =0;i<rd_obs.size();i++){
		try{
		  if(rd_obs[i]>prior_tree_matrix_temp.nrow()){
			throw std::range_error("can't add col because rd row index is out of range");
		  }
		}catch(...){
		  ::Rf_error("there is a problem adding rd col to mat");
		}    
		M(rd_obs[i],ncol)=grow_node+2;
	}
	return(wrap(M));
} 

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_tree(int grow_node,NumericMatrix prior_tree_table_temp,double left_daughter){
	int nrow=prior_tree_table_temp.nrow();
	arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
	M(grow_node-1,5)=0;
	M(grow_node-1,6)=0;
	M.insert_rows(nrow,2);
	M(grow_node-1,0)=left_daughter;
	M(grow_node-1,1)=left_daughter+1;
	M(left_daughter-1,4)=-1;
	M(left_daughter,4)=-1;

	NumericMatrix s=as<NumericMatrix>(wrap(M));
	IntegerVector rname=seq_len(s.nrow());

	List dimnms = // two vec. with static names
	List::create(rname,
			   CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
	// and assign it
	s.attr("dimnames") = dimnms;

	return(s);
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_mat(double d,NumericMatrix prior_tree_matrix_temp,double left_daughter,NumericVector ld_obs,NumericVector rd_obs){
	int ncol_mat=prior_tree_matrix_temp.ncol();
	arma::mat N=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
	arma::vec colmat=N.col(d);
	NumericVector colmat2=wrap(colmat);

	if(d+1==ncol_mat){
		//update prior_tree_matrix
		//insert extra column for the new split node

		N.insert_cols(ncol_mat,1); 
		colmat2[ld_obs]=left_daughter;
		colmat2[rd_obs]=left_daughter+1;
		colmat=Rcpp::as<arma::vec>(colmat2);
		N.col(d+1)=colmat;

	}else{
		colmat2[ld_obs]=left_daughter;
		colmat2[rd_obs]=left_daughter+1;
		colmat=Rcpp::as<arma::vec>(colmat2);
		N.col(d)=colmat;  
	}

	return(wrap(N));
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector remove_zero(NumericVector nodes_at_depth){
	arma::vec nodes_at_depth2=Rcpp::as<arma::vec>(nodes_at_depth);
	arma::vec ret=nodes_at_depth2.elem(arma::find(nodes_at_depth2!=0));
	return(wrap(ret));
}

//######################################################################################################################//
  
// [[Rcpp::export]]
IntegerVector order_intvec_(IntegerVector x) {
	IntegerVector sorted = clone(x).sort();
	std::reverse(sorted.begin(), sorted.end());

	return match(sorted, x);
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector get_gnp(NumericVector nodes_at_depth,int grow_node){
	arma::uvec grow_node_pos=arma::find(as<arma::vec>(nodes_at_depth)==grow_node);

	return(wrap(grow_node_pos));  
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_nodes(NumericMatrix tree_table){
	arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
	arma::vec colmat=arma_tree.col(4);
	arma::uvec term_nodes=arma::find(colmat==-1);
	term_nodes=term_nodes+1;
	
	return(wrap(term_nodes));
} 

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_term_obs(NumericMatrix tree_matrix_temp,double terminal_node){
	arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
	arma::uvec term_obs;

	for(int j=0;j<tree_matrix_temp.ncol();j++){
		arma::vec colmat=arma_tree_mat.col(j);
		term_obs=arma::find(colmat==terminal_node);
		if(term_obs.size()>0){
			break;
		}
	}

	return(term_obs);
}

//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double likelihood_function(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
	double tree_log_lik;
	NumericVector terminal_nodes=find_term_nodes(treetable_temp);
	double b=terminal_nodes.size();
	IntegerVector n(b);
	double term1=0;
	double term2=0;
	arma::uvec term_obs;
	int ni;
	arma::vec y_temp2=as<arma::vec>(y_temp);
	for(int i=0;i< b;i++){
		//number of observations in terminal nodes
		term_obs=find_term_obs(obs_to_nodes_temp,terminal_nodes[i]);
		arma::vec y_k=y_temp2.elem(term_obs);
		ni=term_obs.size();
		n[i]=ni;
		double ybar=0;
		if(y_k.size()!=0){
		  ybar=mean(y_k);
		}
		term1+=log(ni+a);
		arma::vec y_k_sq=pow(y_k,2);
		double sum_yksq=sum(y_k_sq);
		double b2=pow(ni*ybar +a*mu,2)/(ni+a);
		term2+=(sum_yksq+a*pow(mu,2)-b2+nu*lambda);
	}
	tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu)/2)*log(term2);
	for(int i=0;i<b;i++){
		if(n[i]<=5){
			tree_log_lik=tree_log_lik-(100000);
		}
		else{
			tree_log_lik=tree_log_lik;
		}
	}

	return(tree_log_lik);
}  
//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
  
arma::uvec find_internal_nodes(NumericMatrix treetable){

	arma::mat arma_tree(treetable.begin(),treetable.nrow(), treetable.ncol(), false); 
	arma::vec colmat=arma_tree.col(4);
	arma::uvec term_nodes=arma::find(colmat==1);
	term_nodes=term_nodes+1;

	return(term_nodes);
}
  
//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double find_prev_nonterm(arma::uvec find_nonterm,NumericVector prev){
  
	double ret=0;
	int z=prev.size();
	for(int j=0;j<z;j++){
		arma::uvec term_equal = arma::find(find_nonterm==prev[j]);    
		ret+=term_equal.size(); 
	}  
  
	return(ret);
}
  
//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_nodes_to_update(arma::uvec all_ld,double left_daughter){
	arma::uvec gr_ld=arma::find(all_ld>=left_daughter);  
	return(gr_ld);
}
  
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_tree_to_middle(NumericVector node_to_update,NumericMatrix prior_tree_table_temp,int grow_node,double left_daughter){
	for(int i=0;i<node_to_update.size();i++){
		if(prior_tree_table_temp(node_to_update[i],0) && prior_tree_table_temp(node_to_update[i],1)!=0){
			prior_tree_table_temp(node_to_update[i],0)+=2;
			prior_tree_table_temp(node_to_update[i],1)+=2;
		}
	}
  
	prior_tree_table_temp(grow_node-1,5)=0;
	prior_tree_table_temp(grow_node-1,6)=0;

	arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
	M.insert_rows(left_daughter-1,2);
	M(left_daughter-1,4)=-1;
	M(left_daughter,4)=-1;      

	M(grow_node-1,0)=left_daughter;
	M(grow_node-1,1)=left_daughter+1;
	NumericMatrix t=as<NumericMatrix>(wrap(M));
	IntegerVector rname=seq_len(t.nrow());
  
	List dimnms = // two vec. with static names
	List::create(rname,
	CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
	// and assign it
	t.attr("dimnames") = dimnms;
  
	return(t);
}
  
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	NumericMatrix update_grow_obs(NumericMatrix prior_tree_matrix_temp,double grow_node,double left_daughter,double d,NumericVector ld_obs,NumericVector rd_obs){
		arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);
		arma::vec ptm2=prior_tree_matrix_temp2.col(d);
		NumericVector ptm=wrap(ptm2);
		ptm[ld_obs]=left_daughter;
		ptm[rd_obs]=left_daughter+1;
		prior_tree_matrix_temp(_,d)=ptm;  

		return(prior_tree_matrix_temp);
	}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
  
NumericMatrix find_obs_to_update_grow(NumericMatrix prior_tree_matrix_temp,double left_daughter,double d,NumericVector ld_obs,NumericVector rd_obs){
  
	int rows=prior_tree_matrix_temp.nrow();
	int cols=prior_tree_matrix_temp.ncol();
	int elements=rows*cols;
	std::vector<double> rows_obs(elements);
	std::vector<double> cols_obs(elements);
	int count=0;
	for(int i=0;i<prior_tree_matrix_temp.nrow();i++){
		for(int j=0;j<prior_tree_matrix_temp.ncol();j++){
			if(prior_tree_matrix_temp(i,j)>=left_daughter){
				rows_obs[count]=i;
				cols_obs[count]=j;
				count++;
			}
		}
	}
	rows_obs.resize(count);
	cols_obs.resize(count);
  
	if(rows_obs.size()!=0){
		for(int k=0;k< count;k++){
			if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])<left_daughter){
				prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=prior_tree_matrix_temp(rows_obs[k],cols_obs[k]);
			}else if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])==0){
				prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=0;
			}else{   
				int temp=prior_tree_matrix_temp(rows_obs[k],cols_obs[k])+2;
				prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=temp;
			}
		}
	}
	//update prior_tree_matrix
	//insert extra column for the new split node 
	if(prior_tree_matrix_temp.ncol()>d+1){
		arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
		M.insert_cols(prior_tree_matrix_temp.ncol(),1);  
		NumericMatrix prior_tree_matrix=as<NumericMatrix>(wrap(M));
	}
  
	arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);
	arma::vec ptm2=prior_tree_matrix_temp2.col(d+1);
	NumericVector ptm=wrap(ptm2);
	ptm[ld_obs]=left_daughter;
	ptm[rd_obs]=left_daughter+1;
	prior_tree_matrix_temp(_,d+1)=ptm;  

	return(prior_tree_matrix_temp);
}
//######################################################################################################################// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_subset(arma::mat& xmat,NumericVector grow_obs){
	int p = xmat.n_cols;
	arma::mat B(grow_obs.size(),p);

	for(int i=0;i<grow_obs.size();i++){
		B.row(i)=xmat.row(i);
	}

	return(B);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_daughter_obs(arma::mat& xmat,NumericVector obs_to_update,int split_var,double split_point){

	arma::uvec ld_obs_ind(obs_to_update.size());
	arma::uvec rd_obs_ind(obs_to_update.size());

	arma::colvec sv_col;
	List daughter_obs(2);

	sv_col=xmat.col(split_var-1);

	arma::uvec obs_to_update_arma=as<arma::uvec>(obs_to_update);
	ld_obs_ind = arma::find(sv_col.elem(obs_to_update_arma)<=split_point);
	rd_obs_ind = arma::find(sv_col.elem(obs_to_update_arma)>split_point);

	NumericVector ld_ind2(as<NumericVector>(wrap(ld_obs_ind)));
	NumericVector rd_ind2(as<NumericVector>(wrap(rd_obs_ind)));

	NumericVector ld_obs=obs_to_update[ld_ind2];
	NumericVector rd_obs=obs_to_update[rd_ind2];

	daughter_obs[0]=ld_obs;
	daughter_obs[1]=rd_obs;

	return(daughter_obs);

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_cols(NumericMatrix tree_matrix_temp,int terminal_node){

	arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);
	int count=0;
	std::vector<double> term_cols(tree_matrix_temp.ncol());

	for(int j=0;j<tree_matrix_temp.ncol();j++){

		arma::vec tempcol=tree_matrix_temp2.col(j);
		arma::uvec term_nodes=find(tempcol==terminal_node);

		if(term_nodes.size()>0){  
			term_cols[count]=j;
			count++;
		}
	}

	term_cols.resize(count);

	return(wrap(term_cols));

}
  
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_grow_obs(arma::mat& xmat,NumericVector grow_obs,int split_var){

	arma::vec sv_col=xmat.col(split_var-1);
	arma::uvec grow_obs2(as<arma::uvec>(grow_obs));
	arma::vec get_min=sv_col.elem(grow_obs2);

	return(wrap(get_min));
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
  
	List grow_tree(arma::mat& xmat,NumericVector y,NumericMatrix prior_tree_matrix,int grow_node,NumericMatrix prior_tree_table,int splitvar,double splitpoint,
	NumericVector terminal_nodes,NumericVector grow_obs,double d,NumericVector get_min,arma::mat& data_curr_node)
	{

	NumericMatrix prior_tree_matrix_temp=clone(prior_tree_matrix);
	NumericMatrix prior_tree_table_temp=clone(prior_tree_table);
	double yy=xmat.n_cols;
	IntegerVector xx=seq_len(yy);
	prior_tree_table_temp(grow_node-1,3)=splitpoint;
	prior_tree_table_temp(grow_node-1,2)=splitvar;
	prior_tree_table_temp(grow_node-1,4)=1;
	//get data subset for left and right daughter nodes
	List daughter_obs=get_daughter_obs(xmat,grow_obs,splitvar,splitpoint);
	NumericVector ld_obs=daughter_obs[0];
	NumericVector rd_obs=daughter_obs[1];

	if(prior_tree_table_temp.nrow()==grow_node){
		prior_tree_table_temp=add_rows(prior_tree_table_temp,grow_node);
		prior_tree_matrix_temp=addcol(prior_tree_matrix_temp,grow_node,ld_obs,rd_obs);  
	}else{
		//if grow node is in the middle of the tree
		NumericVector nodes_d;
		nodes_d=prior_tree_matrix_temp(_,d);
		NumericVector nodes_at_depth=sort_unique(nodes_d);
		NumericVector nodes_at_depth1=remove_zero(nodes_at_depth);
		NumericVector gn_pos=get_gnp(nodes_at_depth1, grow_node);
		arma::uvec prev_uvec= arma::find(as<arma::vec>(nodes_at_depth1)<grow_node);
		arma::vec nontermvec=Rcpp::as<arma::vec>(nodes_at_depth1);
		nontermvec=nontermvec.elem(prev_uvec);
		NumericVector prev= as<NumericVector>(wrap(nontermvec));
		double prev_nonterm=0;
		if(prev.size()!=0){
			arma::uvec find_nonterm=find_internal_nodes(prior_tree_table);
			//should only find internal nodes at the current depth
			prev_nonterm=find_prev_nonterm(find_nonterm,prev);
		}
		double left_daughter=grow_node +2*(prev_nonterm)+(nodes_at_depth1.size()-gn_pos[0]);
		NumericVector ptt=prior_tree_table(_,1);
		arma::uvec node_to_update=find_nodes_to_update(as<arma::uvec>(ptt),left_daughter);
		//increase the node number of nodes after the grow node by two (because 2 daughter nodes are now appended to grow node)
		//do this for all observations except those that already belong to a terminal node (a value of 0)
		if(node_to_update.size()==0){
			if(prior_tree_matrix_temp.ncol()>d+1){
				//  std::cout<<"add to the end of the table but tree depth isn't increased "<<"\n";
				prior_tree_table_temp=set_daughter_to_end_tree(grow_node,prior_tree_table_temp,left_daughter);
				prior_tree_matrix_temp=update_grow_obs(prior_tree_matrix_temp,grow_node,left_daughter,d+1,ld_obs,rd_obs);
			}else{
				//if the daughter node number already exists in the tree and existing node numbers have to be updated
				//daughter nodes need to be added to the end of the table not in the center of it
				prior_tree_table_temp=set_daughter_to_end_tree(grow_node,prior_tree_table_temp,left_daughter);
				prior_tree_matrix_temp=set_daughter_to_end_mat(d,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);
			}
		}else{
			//if the daughter node number already exists in the tree and existing node numbers have to be updated
			prior_tree_table_temp=set_tree_to_middle(wrap(node_to_update),prior_tree_table_temp,grow_node,left_daughter);
			prior_tree_matrix_temp=find_obs_to_update_grow(prior_tree_matrix_temp,left_daughter,d,ld_obs,rd_obs);    
		}
	}

	List ret(2);
	ret[0]=prior_tree_table_temp;
	ret[1]=prior_tree_matrix_temp;

	return(ret);
}
  
  //######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,NumericMatrix tree_matrix_temp,double term_cols){

	arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);
	arma::vec arma_col=tree_matrix_temp2.col(term_cols+1);
	NumericVector col(as<NumericVector>(wrap(arma_col)));
	col[ld_obs]=left_daughter;
	col[rd_obs]=right_daughter;  
	tree_matrix_temp(_,term_cols+1)=col;

	return(tree_matrix_temp);
}
  
//######################################################################################################################//
  
// [[Rcpp::export]]

IntegerVector order_(NumericVector x) {
	NumericVector sorted = clone(x).sort();
	std::reverse(sorted.begin(), sorted.end());

	return match(sorted, x);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double get_tree_prior(NumericMatrix tree_table,NumericMatrix tree_matrix,double alpha,double beta){
	double propsplit=1;
	IntegerVector d;
	int col=tree_matrix.ncol();
	std::vector<int> int_nodes_index(100*col);
	int index_count=0;  
	arma::uvec internal_nodes_prop=find_internal_nodes(tree_table);
	arma::mat tree_matrix2(tree_matrix.begin(),tree_matrix.nrow(),tree_matrix.ncol(),false);
	int count=internal_nodes_prop.size();

	for(int k=0;k<count;k++){ 
		for(int j=0;j<tree_matrix.ncol();j++){
			arma::vec armacol=tree_matrix2.col(j);
			arma::uvec found=find(armacol==internal_nodes_prop[k]);      
			if(found.size()>0){        
				int_nodes_index[index_count]=j+1;
				index_count++;
				break;
			}        
		}
		int_nodes_index.resize(index_count);
		if(int_nodes_index.size()!=0){      
			d=unique(as<IntegerVector>(wrap(int_nodes_index)));
			double d1=d[0];
			propsplit*=alpha*pow((d1+1),-beta) ;  
		}
		std::vector<int> temp(col);
		int_nodes_index=temp;
		index_count=0;
	} 

	return(propsplit);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_tree(double start_mean,double start_sd){

	NumericMatrix treemat(1,7);
	double rand=R::rnorm(start_mean,start_sd);
	NumericVector testrow = NumericVector::create(0,0,0,0,-1,rand,0);
	for(int k=0;k<1;k++){
		for(int j=0;j<7;j++){
			treemat(k,j)=testrow[j];
		}
	}
	List dimnms = // two vec. with static names
	List::create(CharacterVector::create("1"),
	CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
	// and assign it
	treemat.attr("dimnames") = dimnms;

	return(treemat);
}
//######################################################################################################################//
  
// [[Rcpp::export]]
NumericMatrix start_matrix(int n){
	NumericMatrix mat(n,1);
	std::fill(mat.begin(), mat.end(), 1);
	return(mat);
}
  
//######################################################################################################################//
  
// [[Rcpp::export]]
List evaluate_model_occams_window(NumericVector tree_lik,double lowest_BIC,int c,List tree_list,List tree_mat_list,IntegerVector tree_parent){
	IntegerVector sorted_lik_index=order_(tree_lik);
	std::vector<double> to_be_removed(tree_lik.size());
	int s=0;

	// check tree is in Occam's window if it isn't then delete it from list                    
	while((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC) > c){
		if(s==(tree_lik.size()-1)){
			break;
		}
		//delete tree from tree list
		if((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC)>c){
			//set indicator for index of trees to be removed 
			to_be_removed[s]=sorted_lik_index[s]-1;
			s+=1;
		}
	}

	to_be_removed.resize(s);
	IntegerVector remove_order_index(to_be_removed.size());
	//delete elements from the higest index down 
	remove_order_index=order_(wrap(to_be_removed));

	for(int j=0;j<s;j++){    
		tree_list.erase(to_be_removed[remove_order_index[j]-1]);
		tree_mat_list.erase(to_be_removed[remove_order_index[j]-1]);
		tree_lik.erase(to_be_removed[remove_order_index[j]-1]);
		tree_parent.erase(to_be_removed[remove_order_index[j]-1]);
	}
	List ret(4);
	ret(0)=tree_lik;
	ret(1)=tree_list;
	ret(2)=tree_mat_list;
	ret(3)=tree_parent;

	return(ret);  
}
  
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_testdata_term_obs(NumericMatrix test_data,NumericMatrix tree_data,NumericVector term_node_means) {
	//Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
	//for each tree accepted in Occam's Window.

	arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
	arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
	NumericVector terminal_nodes=find_term_nodes(tree_data);
	arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
	NumericVector tree_predictions;
	//for each internal node find the observations that belong to the terminal nodes
	NumericVector predictions(test_data.nrow());

	for(int i=0;i<terminal_nodes.size();i++){
		arma::mat subdata=testd;
		int curr_term=terminal_nodes[i];
		int row_index;
		int term_node=terminal_nodes[i];
		if(curr_term % 2==0){
			row_index=terminal_nodes[i];
		}else{
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
			node_split_mat.insert_rows(0,1);
			node_split_mat(0,0)=tree_data(parent_node[0],2);
			node_split_mat(0,1)=tree_data(parent_node[0],3);
			node_split_mat(0,2)=rd;
			row_index=parent_node[0] +1;
			term_node=parent_node[0]+1;
		}  
		//fill in the predicted value for tree
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
	} 

	return(predictions);
}
//######################################################################################################################//
  
// [[Rcpp::export]]
List resize(const List& x, int n ){
	List y(n) ;
	for( int i=0; i<n; i++) y[i] = x[i] ;

	return y ;
}
//######################################################################################################################//
  
// [[Rcpp::export]]
List resize_bigger( const List& x, int n ){
	int oldsize = x.size() ;
	List y(n) ;
	for( int i=0; i<oldsize; i++) y[i] = x[i] ;
	return y ;
}
//######################################################################################################################//
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split(NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,double a,double mu,double nu,double lambda,int c,double lowest_BIC,int parent,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize){
	//this function will search through all predictive split points and return those within Occam's Window.
	//std::cout<<"in first round"<<"\n"; 
	int split_var;
	NumericMatrix treetable_c=treetable;
	NumericMatrix treemat_c=tree_mat;

	NumericVector terminal_nodes=find_term_nodes(treetable_c);
	IntegerVector change_node1;
	int list_size=1000;
	std::vector<double> tree_lik(list_size);
	List proposal_tree;
	List ret(9);
	bool no_tree_err=0;
	List likeliest_tree;
	List tree_list(list_size);
	List tree_mat_list(list_size);
	int count=0;
	std::vector<int> tree_parent(list_size);
	int best_sv;
	double best_sp;
	double tree_prior;
	List changetree;
	double BIC;
	int p;
	List eval_model;
	NumericVector int_nodes;
	arma::colvec curr_col=data.col(0);
	arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
	NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
	arma::mat data_curr_node=get_subset(data,wrap(grow_obs));
	double d=d1[0];
	NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
	double lik;

	for(int l=0;l<terminal_nodes.size();l++){
		//loop over each terminal node
		grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
		//depth of tree at current terminal node
		d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
		data_curr_node=get_subset(data,wrap(grow_obs));
		d=d1[0];
		int w=cp_mat.nrow();
		if(data_curr_node.n_rows<=5){
			throw std::range_error("not enough obs in node to grow any further");
			//continue;
		}
		for(int k=0;k<w;k++){
			split_var=cp_mat(k,0)+1;
			arma::colvec curr_cols=data.col(split_var-1);
			get_min=get_grow_obs(data,wrap(grow_obs),split_var);

			if(get_min.size()<=5){
				throw std::range_error("obs in this terminal node are too small");
			}

			double split_point=cp_mat(k,1);
			arma::vec curr_cols2=data_curr_node.col(split_var-1);

			arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
			arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));

			if(ld_prop.size()<=5 || rd_prop.size()<=5){
				continue;
			}
			proposal_tree=grow_tree(data,resids,treemat_c,terminal_nodes[l],treetable_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node);
			NumericMatrix test =proposal_tree[0];
			NumericMatrix test1 =proposal_tree[1];

			if(test1.ncol()==3){
				NumericVector u1=unique(test1(_,0));
				NumericVector u2=unique(test1(_,1));
				NumericVector u3=unique(test1(_,2));
			}
			lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);  
			tree_prior=get_tree_prior(proposal_tree[0],proposal_tree[1],alpha,beta);
			int_nodes=find_term_nodes(proposal_tree[0]);
			p=int_nodes.size();
			BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows);  
			if(BIC<lowest_BIC){
				lowest_BIC=BIC;
				best_sv=split_var;
				best_sp=split_point;
				likeliest_tree=proposal_tree;
				tree_list[count]=proposal_tree[0];
				tree_mat_list[count]=proposal_tree[1];
				tree_lik[count]=BIC;
				tree_parent[count]=parent;
				count++;
				if(count==(tree_list.size()-1)){
					list_size=list_size*2;
					tree_list=resize_bigger(tree_list,list_size);
					tree_mat_list=resize_bigger(tree_mat_list,list_size);
					tree_lik.resize(list_size);
					tree_parent.resize(list_size);
				}
			}else{
				if((BIC)-(lowest_BIC)<=c){
					if(is<NumericMatrix>(proposal_tree[0])){
						//std::cout<<"its a matrix "<<"\n";
					}else{
						throw std::range_error("proposal tree not a matrix");
					}
					tree_list[count]=proposal_tree[0];
					tree_mat_list[count]=proposal_tree[1];
					tree_lik[count]=BIC;
					tree_parent[count]=parent;
					count++;
					if(count==(tree_list.size()-1)){
						list_size=list_size*2;
						tree_list=resize_bigger(tree_list,list_size);
						tree_mat_list=resize_bigger(tree_mat_list,list_size);
						tree_lik.resize(list_size);
						tree_parent.resize(list_size);						  
					}
				}
			}
		}  
	}
	tree_list=resize(tree_list,count);
	tree_mat_list=resize(tree_mat_list,count);
	tree_lik.resize(count);
	tree_parent.resize(count);
	if(count>0){
		eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
		NumericVector testlik =eval_model[0];
		List testtree =eval_model[1];    
		List testmat =eval_model[2]; 
		IntegerVector testpar =eval_model[3];

		if(testlik.size()>0){
			//check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
			if(testlik.size()>maxOWsize){
				IntegerVector owindices=order_(testlik);
				owindices=owindices-1;
				//get the top maxOWsize indices to keep in OW
				NumericVector temp_olik(maxOWsize);
				List temp_otrees(maxOWsize);
				List temp_omat(maxOWsize);
				IntegerVector temp_oparent(maxOWsize);
				for(int t=0;t<maxOWsize;t++){
					temp_olik[t]=testlik[owindices[t]];
					temp_otrees[t]=testtree[owindices[t]];
					temp_omat[t]=testmat[owindices[t]];
					temp_oparent[t]=testpar[owindices[t]];
				}
				testlik=temp_olik;
				testtree=temp_otrees;
				testmat=temp_omat;
				testpar=temp_oparent;
			}
			ret[0]=lowest_BIC;
			ret[1]=best_sv;
			ret[2]=best_sp;
			ret[3]=likeliest_tree;
			ret[4]=testtree;
			ret[5]=testlik;
			ret[6]=testmat;
			ret[7]=testpar;
			ret[8]=no_tree_err;

			return (ret);
		}else{
			//if no trees are found within Occam's window function will return an error to main
			no_tree_err=1;
			List gr(1);
			gr[0]=no_tree_err;
			return(gr);
		}
	}else{
		no_tree_err=1;
		List gr(1);
		gr[0]=no_tree_err;
		return(gr);
	}
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector update_mean_var(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector resids,double a){
	List update_params(1);
	NumericVector terminal_nodes;
	arma::uvec term_obs;
	terminal_nodes= find_term_nodes(tree_table);
	NumericVector Tj(terminal_nodes.size());
	NumericVector new_mean(terminal_nodes.size());
	arma::vec armaresids=as<arma::vec>(resids);

	for(int k=0;k< terminal_nodes.size();k++){    
		term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);
		//get the number of observations in node k
		Tj[k]=term_obs.size();
		NumericVector  get_mean(term_obs.size());
		for(int i=0;i<Tj[k];i++){
			get_mean[i]=resids[term_obs[i]];
		} 
		double sum_resids=std::accumulate(get_mean.begin(),get_mean.end(),0.0);
		new_mean[k]=sum_resids/(Tj[k]+a);
		arma::uvec temp;
		term_obs=temp;
	}  

	return(new_mean);
}
//######################################################################################################################//

// [[Rcpp::export]]
List update_predictions(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector new_mean,int n){

	List updated_preds(2);
	NumericVector new_preds(n);
	NumericVector terminal_nodes;
	arma::uvec term_obs;
	terminal_nodes=find_term_nodes(tree_table);

	for(int k=0;k<terminal_nodes.size();k++){
		term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);        
		//update the terminal node mean of the selected tree nodes:
		tree_table(terminal_nodes[k]-1,5)= new_mean[k];
		IntegerVector term_obs2=wrap(term_obs);
		new_preds[term_obs2]=new_mean[k];
	}
	updated_preds[0]=tree_table;
	updated_preds[1]=new_preds;

	return(updated_preds);
}
//######################################################################################################################//

using namespace Rcpp;
using namespace std;

const double flagval = __DBL_MIN__; 
inline double flag(double a, bool b) { return b ? a : flagval; }

// [[Rcpp::export]]
NumericVector subsetter(NumericVector a, LogicalVector b) {
	NumericVector a1=clone(a);
	transform(a1.begin(), a1.end(), b.begin(), a1.begin(), flag);
	NumericVector res = NumericVector(sum(b));  
	remove_copy(a1.begin(), a1.end(), res.begin(), flagval);

	return res;    
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector order_inc_(NumericVector x) {
	NumericVector sorted = clone(x).sort();
	return match(sorted, x);
}
    
//######################################################################################################################//

// [[Rcpp::export]]
List min_which2(NumericVector array,int n,double minout,int whichout){
	// Function to find minimum of an array with n elements that is put in min
	minout=array[0];
	whichout=0;

	for(int i=1;i<n;i++){
		if(array[i]< minout){
			minout= array[i];
			whichout=i;
		}
	}
	List ret(2);
	ret[0]=minout;
	ret[1]=whichout;

	return(ret);
}
//######################################################################################################################//

#include <Rmath.h>
// [[Rcpp::export]]
double mll_meanvar2(double x, double x2, int n){
	double sigsq=(x2-((x*x)/n))/n;
	if(sigsq<=0){sigsq=0.00000000001;}

	return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector PELT_meanvar_norm2(NumericVector resp,double pen)
// 0 by default, nonzero indicates error in code 
{
	int n=resp.size();
	NumericVector y2=cumsum(pow(resp,2));
	y2.push_front(0);
	NumericVector y=cumsum(resp);
	y.push_front(0);
	IntegerVector cptsout(n,0);    
	IntegerVector lastchangecpts(2*(n+1));
	NumericVector lastchangelike(n+1);
	IntegerVector checklist(n+1);  
	int nchecklist;
	double minout;
	NumericVector tmplike(n+1);
	IntegerVector tmpt(n+1);
	int tstar,i,whichout,nchecktmp;
	double mll_meanvar();
	void min_which();
	lastchangelike[0]= -pen;
	lastchangecpts[0]=0; lastchangecpts[n]=0;
	double x=y[1];
	double x2=y2[1];
	lastchangelike[1]=mll_meanvar2(x,x2,1);
	lastchangecpts[1]=0; lastchangecpts[n+1]=1;
	lastchangelike[2]=mll_meanvar2(y[2],y2[2],2);
	lastchangecpts[2]=0; lastchangecpts[n+2]=2;
	lastchangelike[3]=mll_meanvar2(y[3],y2[3],3);
	lastchangecpts[3]=0; lastchangecpts[n+3]=3;

	minout=lastchangelike[checklist[0]] + mll_meanvar2(x,x2,0)+pen;
	whichout=0;

	nchecklist=2;
	checklist[0]=0;
	checklist[1]=2;

	for(tstar=4;tstar<(n+1);tstar++){
		R_CheckUserInterrupt(); // checks if user has interrupted the R session and quits if true 

		for(i=0;i<nchecklist;i++){
			tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar2(y[tstar]- y[checklist[i]],y2[tstar]-y2[checklist[i]],tstar-checklist[i])+pen;
		}
		List mw=min_which2(tmplike,nchecklist,minout,whichout); //updates minout and whichout with min and which element 
		NumericVector tempmin=mw[0];
		minout=tempmin[0];
		lastchangelike[tstar]=minout;
		whichout=mw[1];
		lastchangecpts[tstar]=checklist[whichout]; 
		lastchangecpts[n+tstar]=tstar;

		// Update checklist for next iteration, first element is next tau 
		nchecktmp=0;
		for(i=0;i<nchecklist;i++){
			if(tmplike[i]<= (lastchangelike[tstar]+pen)){
				checklist[nchecktmp]=checklist[i];
				nchecktmp+=1;
			}
		}
		checklist[nchecktmp]=tstar-1;  // atleast 2 obs per seg
		nchecktmp+=1;
		nchecklist=nchecktmp;
	} // end taustar

	// put final set of changepoints together
	int ncpts=0;
	int last=n;
	while(last!=0){
		cptsout[ncpts]=lastchangecpts[n+last];
		last=lastchangecpts[last];
		ncpts+=1;
	}

	IntegerVector cptsoutret=cptsout[cptsout>0];
	std::sort(cptsoutret.begin(), cptsoutret.end());

	return(cptsoutret);
}
    
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double SS(arma::vec x, arma::vec y, double split){            
	double meanLeft, meanRight;
	int n = x.n_rows;
	arma::vec meanAll(n);
	arma::vec xLeft = y.elem(arma::find(x <= split));
	arma::vec xRight = y.elem(arma::find(x > split));
	meanLeft = mean(xLeft);
	meanRight = mean(xRight);

	for(int j=0;j<n;j++) {

		if(x[j]<=split){
			meanAll[j] = meanLeft;
		}else{
			meanAll[j] = meanRight;
		}  
	}
	arma::vec test_resids=y-meanAll;
	double tSS=as<double>(wrap(trans(y-meanAll)*(y-meanAll)));

	return tSS;
}
//######################################################################################################################//
// [[Rcpp::export]]

List gridCP(arma::vec x, arma::vec y, int gridSize = 10) {

	NumericVector out(gridSize-2);
	NumericVector cp_strength(gridSize-2);
	arma::vec no_split(gridSize-2);
	double xGrid, gridStep = (max(x)-min(x))/((double)gridSize-1.0);
	xGrid = min(x);
	List currSS(2);

	for(int i=1;i<(gridSize-1);i++) {
		xGrid += gridStep;
		arma::vec ld_size= y.elem(arma::find(x <= xGrid));
		arma::vec rd_size=y.elem(arma::find(x > xGrid));    

		if(ld_size.size()>5 && rd_size.size()>5)
		{
			out[i-1]=xGrid;
			double testSS=SS(x,y,xGrid);
			cp_strength[i-1] = testSS;
			no_split[i-1]=0;
		}else{
			no_split[i-1]=1;
		}
	}
	arma::uvec to_remove=find(no_split ==1);
	IntegerVector remove_order_index=as<IntegerVector>(wrap(order_(as<NumericVector> (wrap(to_remove)))));

	if(to_remove.size()>0){
		for(int k=0;k<remove_order_index.size();k++){
			out.erase(to_remove[remove_order_index[k]-1]);
			cp_strength.erase(to_remove[remove_order_index[k]-1]);
		}
	}
	if(out.size()>0){
		currSS[0]= out;
		currSS[1]= cp_strength;

		return currSS;
	}else{
		List ret(2);
		ret[0]=0;
		ret[1]=0;

		return ret;
	}
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_gridpoint_cpmat(NumericMatrix data,NumericVector resp,int gridsize,int num_cp){
	int err1=0;
	int numvar = data.ncol();
	int numrows=gridsize*numvar-1;
	int numcols=3;
	arma::mat change_points(numrows,numcols);
	int row_count=-1;
	IntegerVector sorted_var;

	for(int i=0;i<data.ncol();i++){
		NumericVector coli=data(_,i);
		arma::vec coli_ar=as<arma::vec>(coli);   
		arma::vec t=Rcpp::as<arma::vec>(resp);
		List ans1 = gridCP(coli,t,gridsize);
		NumericVector ans2=ans1[0];
		NumericVector cp_strength=ans1[1];

		if(ans1.size()!=1 && ans2[0]!=0){
			for(int j=0;j<ans2.size();j++){
				row_count+=1;
				change_points(row_count,0)=i;
				change_points(row_count,1)=ans2[j];
				change_points(row_count,2)=cp_strength[j];
			}   
		}
	}

	if(row_count+1!=(int) change_points.n_rows){
		change_points.shed_rows(row_count+1,change_points.n_rows-1);
	}
	arma::vec te=change_points.col(2);
	NumericVector col_to_order=as<NumericVector>(wrap(te));
	IntegerVector ordered_dev=order_inc_(col_to_order);
	ordered_dev=ordered_dev-1;
	change_points.shed_col(2);
	int cp=change_points.n_rows;
	double cp_prop=(double)num_cp/(double)100;
	int num_cp2=round(cp*(cp_prop));
	num_cp=round(num_cp2);

	if(num_cp==0 && cp!=0){
		num_cp=cp;
	}
	if(cp<num_cp){
		num_cp=change_points.n_rows;
	}
	if(num_cp==0){
		err1=1;
	}
	if(err1==0){
		arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);
		t=t.subvec(0,num_cp-1);
		change_points=change_points.rows(t);
	}
	List ret(2);
	ret[0]=wrap(change_points);
	ret[1]=err1;

	return(ret);
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_pelt_cpmat(NumericMatrix data,NumericVector resp,double pen,int num_cp){
	int err1=0;
	int n=data.nrow();
	arma::mat change_points;
	change_points.set_size(10000000,3);
	int row_count=-1;
	IntegerVector sorted_var;

	for(int i=0;i<data.ncol();i++){
		NumericVector coli=data(_,i);
		arma::vec coli_ar=as<arma::vec>(coli);
		sorted_var=order_inc_(coli);
		sorted_var=sorted_var-1;
		IntegerVector ans = PELT_meanvar_norm2(resp[sorted_var],pen);
		ans=ans-1;
		arma::vec resp_ar=as<arma::vec>(resp);
		NumericVector allMean(n);

		if(ans.size()!=1 && ans[0]!=n){
			double meanRight;
			double meanLeft;      
			IntegerVector pelt_sp=sorted_var[ans];
			NumericVector split_values=coli[pelt_sp];
			NumericVector unique_sp=unique(split_values);
			NumericVector all_sp=coli[pelt_sp];

			if(unique_sp.size()<all_sp.size()){
				IntegerVector index=match(unique_sp,all_sp);
				arma::vec indexarma=as<arma::vec>(index);
				IntegerVector t4=ifelse(is_na(index),1,0);
				arma::vec t4arma=as<arma::vec>(t4);
				arma::uvec t42=find(t4arma==0);
				IntegerVector index2(t42.size());

				arma::vec index3=indexarma.elem(t42);

				index2=wrap(index3);
				index2=index2-1;
				ans=index2;
			}

			for(int j=0;j<ans.size()-1;j++){
				arma::vec xLeft;
				arma::vec xRight;
				double splitpoint;
				if(unique_sp.size()<all_sp.size()){
					splitpoint = all_sp[j];
					xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));
					xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));
				}else{
					splitpoint = coli[sorted_var[ans[j]]];
					xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));
					xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));
				}
				if(xLeft.size()>5 && xRight.size()>5)
				{
					row_count +=1;
					meanLeft = mean(xLeft);
					meanRight = mean(xRight);
					arma::uvec left_ind=find(coli_ar<=splitpoint);
					arma::uvec right_ind=find(coli_ar>splitpoint);
					NumericVector left_ind2=wrap(left_ind);
					NumericVector right_ind2=wrap(right_ind);
					allMean[left_ind2]=meanLeft;
					allMean[right_ind2]=meanRight;
					arma::vec allMean_ar=as<arma::vec>(allMean);

					double cp_strength;
					cp_strength=as<double>(wrap(trans(resp_ar-allMean_ar)*(resp_ar-allMean_ar)));
					change_points(row_count,0)=i;
					change_points(row_count,1)=coli[sorted_var[ans[j]]];
					change_points(row_count,2)=cp_strength;
				}
			}          
		}
	}
	change_points.shed_rows(row_count+1,change_points.n_rows-1);
	arma::vec te=change_points.col(2);
	NumericVector col_to_order=as<NumericVector>(wrap(te));
	IntegerVector ordered_dev=order_inc_(col_to_order);
	ordered_dev=ordered_dev-1;
	change_points.shed_col(2);
	int cp=change_points.n_rows;
	double cp_prop=(double)num_cp/(double)100;
	int num_cp2=round(cp*(cp_prop));
	num_cp=round(num_cp2);
	if(num_cp==0 && cp!=0){
		num_cp=cp;
	}
	if(cp<num_cp){
		num_cp=change_points.n_rows;
	}
	if(num_cp==0){
		err1=1;
	}
	if(err1==0){
		arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);
		t=t.subvec(0,num_cp-1);
		change_points=change_points.rows(t);
	}
	List ret(2);
	ret[0]=wrap(change_points);
	ret[1]=err1;

	return(ret);
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees(arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,int c,double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,
int first_round,IntegerVector parent,List cp_mat_list,IntegerVector err_list,NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize
){
	List eval_model;
	NumericVector lik_list;
	List best_subset;
	int overall_size=1000;
	List overall_trees(overall_size);
	NumericVector overall_lik2;
	IntegerVector overall_parent2;
	List overall_mat(overall_size);
	int overall_count=0;  
	std::vector<int> overall_parent(overall_size);
	std::vector<double> overall_lik(overall_size);
	NumericVector test_preds;

	for(int j=0;j<5;j++){
		int lsize=1000;
		List table_subset_curr_round(lsize);
		std::vector<double> lik_subset_curr_round(lsize);
		List mat_subset_curr_round(lsize);
		std::vector<int> parent_curr_round(lsize);
		int count=0;
		for(int i=0;i<tree_table.size();i++){
			if(first_round==1){
				parent=-1;
				NumericMatrix temp_list=cp_mat_list[0];
				best_subset=get_best_split(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize); 
			}else{
				if(err_list[i]==0){
					NumericMatrix test_tree=tree_table[i];
					NumericMatrix test_treemat=tree_mat[i];
					NumericMatrix test_cpmat= cp_mat_list[parent[i]];
					best_subset=get_best_split(resids(_,parent[i]),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[i],cp_mat_list[parent[i]],alpha,beta,maxOWsize);    
					}else if(err_list[i]==1){
					continue;
					}else{
					List ret_list(6);
					ret_list[0]=9999;
					ret_list[1]=err_list[i];
					ret_list[2]=i;
					ret_list[3]=j;
					ret_list[4]=tree_table;
					ret_list[5]=err_list;
					return(ret_list);
					throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
				}
			}    
			if(best_subset.size()==1){
				continue;
			}
			List temp_trees=best_subset[4];
			List temp_mat=best_subset[6];
			lik_list=best_subset[5];
			IntegerVector temp_parent=best_subset[7];
			if(temp_parent.size()!= temp_trees.size()){
				throw std::range_error("there should be a parent for each tree!!!");
			}
			if(lik_list.size()==0){
				throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
			}

			if(min(lik_list)<lowest_BIC){
				lowest_BIC=min(lik_list);
			}
			for(int k=0;k<temp_trees.size();k++){
				table_subset_curr_round[count]=temp_trees[k];
				lik_subset_curr_round[count]=lik_list[k];
				mat_subset_curr_round[count]=temp_mat[k];
				parent_curr_round[count]=temp_parent[k];
				count++;

				if(count==(lsize-1)){
					lsize=lsize*2;
					table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
					mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
					lik_subset_curr_round.resize(lsize);
					parent_curr_round.resize(lsize);
				}
			}
		}
		table_subset_curr_round=resize(table_subset_curr_round,count);
		mat_subset_curr_round=resize(mat_subset_curr_round,count);
		lik_subset_curr_round.resize(count);
		parent_curr_round.resize(count);

		if(table_subset_curr_round.size()==0){
			break;
		}
		for(int k=0;k<table_subset_curr_round.size();k++){
			overall_trees[overall_count]=table_subset_curr_round[k];
			overall_lik[overall_count]=lik_subset_curr_round[k];
			overall_mat[overall_count]=mat_subset_curr_round[k];
			overall_parent[overall_count]=parent_curr_round[k];
			overall_count++;

			if(overall_count==(overall_size-1)){
				overall_size=overall_size*2;
				overall_trees=resize_bigger(overall_trees,overall_size);
				overall_lik.resize(overall_size);
				overall_mat=resize_bigger(overall_mat,overall_size);
				overall_parent.resize(overall_size);
			}
		}
		overall_trees=resize(overall_trees,overall_count);
		overall_lik.resize(overall_count);
		overall_mat=resize(overall_mat,overall_count);
		overall_parent.resize(overall_count);
		eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
		overall_lik2=eval_model[0];
		overall_trees=eval_model[1];
		overall_mat=eval_model[2];
		overall_count=overall_trees.size();
		overall_parent2=eval_model[3];
		//add in check to see if OW accepted more than the top maxOW models...
		if(overall_lik2.size()>maxOWsize){
			//find the maxOWsize best models and continue with those!
			IntegerVector owindices=order_(overall_lik2);
			owindices=owindices-1;
			//get the top maxOWsize indices to keep in OW
			NumericVector temp_olik(maxOWsize);
			List temp_otrees(maxOWsize);
			List temp_omat(maxOWsize);
			IntegerVector temp_oparent(maxOWsize);

			//now only select those elements
			for(int t=0;t<maxOWsize;t++){  
				temp_olik[t]=overall_lik2[owindices[t]];
				temp_otrees[t]=overall_trees[owindices[t]];
				temp_omat[t]= overall_mat[owindices[t]];
				temp_oparent[t]=overall_parent2[owindices[t]];
			}

			overall_lik2=temp_olik;
			overall_trees=temp_otrees;
			overall_mat=temp_omat;
			overall_count=overall_trees.size();
			overall_parent2=temp_oparent;
		}

		tree_table=table_subset_curr_round;
		IntegerVector temp1(table_subset_curr_round.size(),1);
		err_list=temp1;
		if(overall_trees.size()<overall_size-1){
			overall_trees=resize_bigger(overall_trees,overall_size);
			overall_mat=resize_bigger(overall_mat,overall_size);
			overall_lik.resize(overall_size);
			overall_parent.resize(overall_size);
		}else{
			overall_size=2*overall_size;
			overall_trees=resize_bigger(overall_trees,overall_size);
			overall_mat=resize_bigger(overall_mat,overall_size); 
			overall_lik.resize(overall_size);
			overall_parent.resize(overall_size);
		}
		tree_mat=mat_subset_curr_round;
		parent=parent_curr_round;

		if(split_rule_node==1){
			NumericVector temp_preds;
			List updated_curr_preds;
			NumericVector new_mean;
			lowest_BIC=min(overall_lik2);
			NumericMatrix curr_resids(resids.nrow(),resids.ncol());

			for(int k=0;k<table_subset_curr_round.size();k++){
				NumericVector terminal_nodes;

				if(parent_curr_round[k]==-1){
					new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
				}else{    
					new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
				}  

				terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
				updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
				NumericVector test_res;

				if(parent_curr_round[k]==-1){
					test_res=resids(_,0); 
				}else{
					test_res=resids(_,parent_curr_round[k]);
				}
				
				NumericVector curr_test_res=updated_curr_preds[1];

				if(parent_curr_round[k]==-1){
					curr_resids(_,0)=test_res-curr_test_res;
				}else{
					curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
				}
			}
			List temp(0);

			cp_mat_list=temp;

			for(int f=0;f<curr_resids.ncol();f++){
				List cp_mat_list1;
			if(gridpoint==0){
				cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
			}else{
				cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
			}

			cp_mat_list.push_back(cp_mat_list1[0]);      
			}

		}
	}
	overall_trees=resize(overall_trees,overall_count);
	overall_mat=resize(overall_mat,overall_count); 
	overall_lik.resize(overall_count);
	overall_parent.resize(overall_count);
	NumericVector temp_preds;
	List updated_preds;
	NumericVector new_mean;  
	NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
	NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
	lowest_BIC=min(overall_lik2);

	for(int k=0;k<overall_trees.size();k++){
		NumericVector terminal_nodes;

		if(overall_parent2[k]==-1){
			new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
		}else{    
			new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
		}  
		terminal_nodes=find_term_nodes(overall_trees[k]);
		updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
		//get the predicted values for the test data.
		if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k],new_mean);
		temp_preds=updated_preds[1];
		overallpreds(_,k)=temp_preds;
		if(is_test_data)overall_test_preds(_,k)=test_preds;
	}
	arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
	arma::colvec predicted_values=sum(M1,1);
	arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
	arma::colvec predicted_test_values=sum(M2,1);
	List ret(8);
	ret[0]=overall_lik2;
	ret[1]=overall_trees;
	ret[2]=overall_mat;
	ret[3]=predicted_values;
	ret[4]=overall_parent2;
	ret[5]=wrap(M1);
	ret[6]=lowest_BIC;
	ret[7]=wrap(M2);
	return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response(double a,double b,double c,double d,NumericVector y){
	NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));

	return(y_scaled);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector get_original(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
	NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);

	return(original_y);
}
//######################################################################################################################//

// [[Rcpp::export]]

List BART_BMA(NumericMatrix data,NumericVector y,double start_mean,double start_sd,double a,double mu,double nu,double lambda,int c,
double sigma_mu,double pen,int num_cp,NumericMatrix test_data,int num_rounds,double alpha,double beta,bool split_rule_node,bool gridpoint,int maxOWsize){
	bool is_test_data=0;
	if(test_data.nrow()>0){
		is_test_data=1;
	}
	if(y.size() !=data.nrow()){
		if(y.size()<data.nrow()){
			throw std::range_error("Response length is smaller than the number of observations in the data"); 
		}else{
			throw std::range_error("Response length is greater than the number of observations in the data"); 
		}
	}
	//check test data has the same number of variables as training data
	if(test_data.nrow()>0 && (data.ncol() != test_data.ncol())){
		throw std::range_error("Test data and training data must have the same number of variables. BART BMA assumes variables are in the same order."); 
	}
	//check value of c is greater than 1!
	if(c<1){
		throw std::range_error("Value of Occam's Window has to be greater than 0."); 
	}
	if(num_cp<0 || num_cp>100){
		throw std::range_error("Value of num_cp should be a value between 1 and 100."); 
	}
	NumericMatrix treetable=start_tree(mu,sigma_mu);
	NumericMatrix treemat=start_matrix(data.nrow());
	//initialize the tree table and matrix
	NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
	//initialize the tree table and matrix
	arma::mat D1(data.begin(), data.nrow(), data.ncol(), false);
	double n=D1.n_rows;
	double lik=likelihood_function(y_scaled,treetable,treemat,a,mu,nu,lambda);
	double tree_prior=get_tree_prior(treetable,treemat,alpha,beta);
	double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);
	List best_subset;  
	List tree_table;
	List tree_mat;
	tree_table.push_back(treetable);
	tree_mat.push_back(treemat);
	List CART_BMA;
	arma::mat r;
	arma::colvec yarma=clone(y_scaled);
	r.insert_cols(0,yarma);
	NumericMatrix resids=wrap(r);
	int first_round;
	List overall_trees(num_rounds);
	List overall_mat;
	List overall_lik;
	NumericMatrix prev_round_preds;  
	NumericVector prev_round_BIC;
	NumericVector prev_round_BIC2;
	arma::mat prev_round_preds2;
	NumericMatrix prev_round_test_preds;
	arma::mat prev_round_test_preds2;
	arma::mat overall_overall_sum_test_preds;
	arma::colvec predicted_test_values;
	List prev_sum_trees;
	List prev_sum_tree_resids;
	List prev_sum_trees_mat;  
	List cp_mat_list;
	int oo_size=300;
	List overall_overall_sum_trees(oo_size);
	List overall_overall_sum_tree_resids(oo_size);
	List overall_overall_sum_trees_mat(oo_size);
	List overall_overall_sum_BIC(oo_size);
	int oo_count=0;
	arma::mat overall_overall_sum_preds;
	IntegerVector prev_par;
	arma::colvec predicted_values;

	for(int j=0;j<num_rounds;j++){
		int overall_size=300;
		List overall_sum_trees(overall_size);
		List overall_sum_trees_mat(overall_size);
		List overall_sum_tree_resids(overall_size);
		int overall_count=0;
		IntegerVector parent;   
		NumericVector curr_round_lik;
		List curr_round_trees;
		List curr_round_mat;
		NumericVector curr_BIC;
		IntegerVector curr_round_parent;
		NumericVector overall_sum_BIC;
		arma::mat overall_sum_preds;    
		arma::mat overall_sum_test_preds;

		if(j==0){
			parent.push_back(0);
			first_round=1;
		}else{
			first_round=0;
		}
		List resids_cp_mat(resids.ncol());
		int resids_count=0;
		std::vector<int> err_list(resids.ncol());

		for(int f=0;f<resids.ncol();f++){
			if(gridpoint==0){
				cp_mat_list=make_pelt_cpmat(data,resids(_,f),pen,num_cp);
			}else{
				cp_mat_list=make_gridpoint_cpmat(data,resids(_,f),pen,num_cp);
			}
			resids_cp_mat[resids_count]=cp_mat_list[0];
			err_list[resids_count]=cp_mat_list[1];
			resids_count++;
		}
		resids_cp_mat=resize(resids_cp_mat,resids_count);
		err_list.resize(resids_count);     
		parent=seq_len(tree_table.size())-1;
		if(is_true(all(as<IntegerVector>(wrap(err_list))==1))){
			if(j==0){
				throw std::range_error("No split points could be found to grow trees");
			}else{
				throw std::range_error("No trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
			}
		}  
		CART_BMA=get_best_trees(D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,first_round,parent,resids_cp_mat,as<IntegerVector>(wrap(err_list)),test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize);

		if(CART_BMA.size()==6 && CART_BMA[0]==9999){
			return(CART_BMA);
		}
		curr_round_lik=CART_BMA[0];
		curr_round_trees=CART_BMA[1];      
		curr_round_mat=CART_BMA[2];
		curr_round_parent=CART_BMA[4];
		NumericMatrix curr_round_preds=CART_BMA[5];
		NumericMatrix curr_round_test_preds=CART_BMA[7];
		curr_BIC=CART_BMA[6];
		if(curr_round_lik.size()==0) {
			break;
		} 

		if(curr_BIC[0]<lowest_BIC){
			lowest_BIC=curr_BIC[0];
		}
		tree_table=List();
		tree_mat=List();
		int lsize=curr_round_lik.size();
		tree_table=List(lsize);
		tree_mat=List(lsize);
		NumericMatrix temp_preds(n,curr_round_lik.size());
		NumericMatrix temp_test_preds(test_data.nrow(),curr_round_lik.size());
		NumericMatrix temp_resids(n,curr_round_lik.size());
		NumericVector temp_parent(curr_round_lik.size());
		NumericVector temp_BIC(curr_round_lik.size());
		List temp_sum_trees(lsize);
		List temp_sum_tree_resids(lsize);
		List temp_sum_trees_mat(lsize);
		int count=0; 
		for(int k=0;k<curr_round_lik.size();k++){
			tree_table[count]=start_tree(mu,sigma_mu);
			tree_mat[count]=start_matrix(n);
			if(j==0){
				temp_preds(_,k)=curr_round_preds(_,k);
				if(is_test_data==1) temp_test_preds(_,k)=curr_round_test_preds(_,k);
				temp_resids(_,k)=y_scaled-temp_preds(_,k);
				temp_parent[k]=-1;
				temp_BIC[k]=curr_round_lik[k];
				temp_sum_trees[count]=curr_round_trees[k];
				temp_sum_trees_mat[count]=curr_round_mat[k];
				temp_sum_tree_resids[count]=resids(_,0);
			}else{
				NumericVector curr_temp_pred=curr_round_preds(_,k) + prev_round_preds(_,curr_round_parent[k]);
				NumericVector curr_temp_test_pred;

				if(is_test_data==1) {
					curr_temp_test_pred=curr_round_test_preds(_,k) + prev_round_test_preds(_,curr_round_parent[k]);
					temp_test_preds(_,k) = curr_temp_test_pred;
				}
				temp_BIC[k]=curr_round_lik[k] + prev_round_BIC[curr_round_parent[k]];
				temp_preds(_,k)=curr_temp_pred;
				temp_resids(_,k)=y_scaled-curr_temp_pred;
				temp_parent[k] = k;
				temp_sum_trees[count]=curr_round_trees[k];
				temp_sum_trees_mat[count]=curr_round_mat[k];
				temp_sum_tree_resids[count]=resids(_,curr_round_parent[k]);        
			}
			count++;     
		}   
		if(curr_round_lik.size()==0){
			throw std::range_error("No trees chosen in last round");
		}
		for(int k=0;k<curr_round_lik.size();k++){
			int size_mat=300;
			List sum_of_trees(size_mat);
			List sum_of_tree_resids(size_mat);
			List sum_of_trees_mat(size_mat);
			int count=0;

			if(curr_round_parent[k]==-1){
			}else{
				if(j==1){
					NumericMatrix other_tree=prev_sum_trees[curr_round_parent[k]];
					NumericVector other_resids=prev_sum_tree_resids[curr_round_parent[k]];
					sum_of_trees[count]= other_tree;  
					sum_of_tree_resids[count]=other_resids;
					NumericMatrix other_mat=prev_sum_trees_mat[curr_round_parent[k]];
					sum_of_trees_mat[count]=other_mat;
					count++;

					if(count==(size_mat-1)){
						size_mat=size_mat*2;
						sum_of_trees=resize_bigger(sum_of_trees,size_mat);
						sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
						sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
					}
				}else{
					List other_tree=prev_sum_trees[curr_round_parent[k]];
					List other_tree_resids=prev_sum_tree_resids[curr_round_parent[k]];  
					List other_mat=prev_sum_trees_mat[curr_round_parent[k]];
					for(int f=0;f<other_tree.size();f++){
						if(is<NumericMatrix>(other_tree[f])){
						}else{
							throw std::range_error(" tree is not a numeric matrix!");
						}
						NumericMatrix treetoadd=other_tree[f];
						if(is<NumericVector>(other_tree_resids[f])){
						}else{
							throw std::range_error("other resids not a numeric matrix!");
						}
						NumericVector residstoadd=other_tree_resids[f];
						if(is<NumericMatrix>(other_mat[f])){
							//  std::cout<<"everything ok"<<"\n";
						}else{
							throw std::range_error(" other mat not a numeric matrix!");
						}
						NumericMatrix mattoadd=other_mat[f];
						sum_of_trees[count]=treetoadd;
						sum_of_tree_resids[count]=residstoadd;
						sum_of_trees_mat[count]=mattoadd;
						count++;

						if(count==(size_mat-1)){
							size_mat=size_mat*2;
							sum_of_trees=resize_bigger(sum_of_trees,size_mat);
							sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
							sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
						}
					}
				}
				sum_of_trees[count]=temp_sum_trees[k];
				sum_of_tree_resids[count]=temp_sum_tree_resids[k];
				sum_of_trees_mat[count]=temp_sum_trees_mat[k];
				count++;

				if(count==(size_mat-1)){
					size_mat=size_mat*2;
					sum_of_trees=resize_bigger(sum_of_trees,size_mat);
					sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
					sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
				}
			}
			sum_of_trees=resize(sum_of_trees,count);
			sum_of_trees_mat=resize(sum_of_trees_mat,count);
			sum_of_tree_resids=resize(sum_of_tree_resids,count);

			if(curr_round_parent[k]!=-1){
				overall_sum_trees[overall_count]=sum_of_trees;
				overall_sum_tree_resids[overall_count]=sum_of_tree_resids;
				overall_sum_trees_mat[overall_count]=sum_of_trees_mat;
				overall_sum_BIC=temp_BIC;
				overall_sum_preds=Rcpp::as<arma::mat>(temp_preds);
				if(is_test_data==1) overall_sum_test_preds=Rcpp::as<arma::mat>(temp_test_preds);
				overall_count++;
				if(overall_count==(overall_size-1)){
					overall_size=overall_size*2;
					overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
					overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
					overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
				}
			}  
		}
		//check if there were any trees from the previous round that didn't have daughter trees grown.
		//create vector to count number of possible parents for previous round
		if(j>0){
			IntegerVector prev_par_no_child=match(prev_par,curr_round_parent);
			if(any(is_na(prev_par_no_child))){
				IntegerVector t4=ifelse(is_na(prev_par_no_child),1,0);
				for(int h=0;h<prev_par_no_child.size();h++){
					if(t4[h]==1){ 
						if(prev_round_BIC2[h]-lowest_BIC<=log(c)){
							SEXP s = prev_sum_trees[h];
							if(is<List>(s)){
								List tree_no_child=prev_sum_trees[h];
								List resids_no_child=prev_sum_tree_resids[h];
								List treemat_no_child=prev_sum_trees_mat[h];
								overall_sum_trees[overall_count]=tree_no_child;
								overall_sum_tree_resids[overall_count]=resids_no_child;
								overall_sum_trees_mat[overall_count]=treemat_no_child;
								overall_count++;

								if(overall_count==(overall_size-1)){
									overall_size=overall_size*2;
									overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
									overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
									overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
								}
								double BIC_to_add=prev_round_BIC2[h];
								overall_sum_BIC.push_back(BIC_to_add);
								overall_sum_preds.insert_cols(overall_sum_preds.n_cols,prev_round_preds2.col(h)); 
								if(is_test_data==1) overall_sum_test_preds.insert_cols(overall_sum_test_preds.n_cols,prev_round_test_preds2.col(h));
							}else{
								NumericMatrix tree_no_child=prev_sum_trees[h];
								NumericVector resids_no_child= prev_sum_tree_resids[h];
								NumericMatrix treemat_no_child=prev_sum_trees_mat[h];
								overall_sum_trees[overall_count]=tree_no_child;
								overall_sum_tree_resids[overall_count]=resids_no_child;
								overall_sum_trees_mat[overall_count]=treemat_no_child;
								overall_count++;

								if(overall_count==(overall_size-1)){
									overall_size=overall_size*2;
									overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
									overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
									overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
								}

								double BIC_to_add=prev_round_BIC2[h];
								overall_sum_BIC.push_back(BIC_to_add);

								overall_sum_preds.insert_cols(overall_sum_preds.n_cols,prev_round_preds2.col(h));                       
								if(is_test_data==1) overall_sum_test_preds.insert_cols(overall_sum_test_preds.n_cols,prev_round_test_preds2.col(h));
							}
						}
					}              
				}
			}          
		}
		//   std::cout<<"Done "<<"\n";
		prev_round_preds=temp_preds;
		if(is_test_data==1) prev_round_test_preds=temp_test_preds;
		prev_round_BIC=temp_BIC;
		prev_round_BIC2=temp_BIC;
		prev_round_preds2=Rcpp::as<arma::mat>(temp_preds);
		if(is_test_data==1) prev_round_test_preds2=Rcpp::as<arma::mat>(temp_test_preds);
		resids=temp_resids;
		parent=temp_parent;    
		overall_sum_trees=resize(overall_sum_trees,overall_count);
		overall_sum_tree_resids=resize(overall_sum_tree_resids,overall_count);
		overall_sum_trees_mat=resize(overall_sum_trees_mat,overall_count);     

		if(first_round==1){  
			prev_sum_trees=temp_sum_trees;
			prev_sum_tree_resids=temp_sum_tree_resids;
			NumericMatrix test=prev_sum_trees[0];
			prev_sum_trees_mat=temp_sum_trees_mat; 
			overall_sum_trees=resize(temp_sum_trees,temp_sum_trees.size());
			overall_sum_trees=temp_sum_trees;
			overall_sum_tree_resids=resize(temp_sum_tree_resids,temp_sum_tree_resids.size());
			overall_sum_tree_resids=temp_sum_tree_resids;
			overall_sum_trees_mat=temp_sum_trees_mat; 
			overall_sum_trees_mat=resize(temp_sum_trees_mat,temp_sum_trees.size());
			overall_sum_BIC=temp_BIC;
			overall_sum_preds= Rcpp::as<arma::mat>(temp_preds);
			if(is_test_data==1) overall_sum_test_preds= Rcpp::as<arma::mat>(temp_test_preds);
		}else{
			prev_sum_trees=overall_sum_trees;
			prev_sum_tree_resids=overall_sum_tree_resids;
			prev_sum_trees_mat=overall_sum_trees_mat;
			prev_round_BIC2=overall_sum_BIC;
			prev_round_preds2=overall_sum_preds;
			if(is_test_data==1) prev_round_test_preds2=overall_sum_test_preds;
		}
		overall_overall_sum_trees[oo_count]=overall_sum_trees;
		overall_overall_sum_tree_resids[oo_count]=overall_sum_tree_resids;
		overall_overall_sum_trees_mat[oo_count]=overall_sum_trees_mat;
		overall_overall_sum_BIC[oo_count]=overall_sum_BIC;
		oo_count ++;

		if(oo_count==(oo_size-1)){
			oo_size=oo_size*2;
			overall_overall_sum_trees=resize_bigger(overall_overall_sum_trees,oo_size);
			overall_overall_sum_tree_resids=resize_bigger(overall_overall_sum_tree_resids,oo_size);
			overall_overall_sum_trees_mat=resize_bigger(overall_overall_sum_trees_mat,oo_size);
			overall_overall_sum_BIC=resize_bigger(overall_overall_sum_BIC,oo_size);
		}    
		overall_overall_sum_preds=overall_sum_preds;

		if(is_test_data==1) overall_overall_sum_test_preds=overall_sum_test_preds;
		overall_trees[j]=curr_round_trees;
		overall_mat.push_back(curr_round_mat);
		overall_lik.push_back(curr_round_lik);
		prev_par=seq_len(overall_sum_trees.size())-1;
	} 
	overall_overall_sum_trees=resize(overall_overall_sum_trees,oo_count);
	overall_overall_sum_tree_resids=resize(overall_overall_sum_tree_resids,oo_count);
	overall_overall_sum_trees_mat=resize(overall_overall_sum_trees_mat,oo_count);
	overall_overall_sum_BIC=resize(overall_overall_sum_BIC,oo_count);
	NumericVector end_BIC=overall_overall_sum_BIC[overall_overall_sum_BIC.size()-1] ;
	NumericMatrix overallpreds(n,end_BIC.size());
	NumericMatrix overall_test_preds(test_data.nrow(),end_BIC.size());
	for(int k=0;k<end_BIC.size();k++){  
		NumericVector temp_preds=prev_round_preds(_,k);
		NumericVector temp_test_preds;
		if(is_test_data==1)temp_test_preds=prev_round_test_preds(_,k);
		NumericVector orig_temp_preds=get_original(min(y),max(y),-0.5,0.5,temp_preds) ;
		NumericVector BICi=-0.5*end_BIC;
		double max_BIC=max(BICi);
		double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
		overallpreds(_,k) = temp_preds*weight;
		if(is_test_data==1)overall_test_preds(_,k) = temp_test_preds*weight;
	}     
	arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
	predicted_values=sum(M1,1);
	arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
	if(is_test_data==1) predicted_test_values=sum(M2,1);
	if(overall_lik.size()==0){
		throw std::range_error("BART-BMA didnt find any suitable model for the data. Maybe limit for Occam's window is too small.");
	}else{
		NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
		NumericVector orig_test_preds;
		if(is_test_data==1){
			orig_test_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_test_values)) ;
		}
		NumericVector minmax(2);
		minmax[0]=min(y);
		minmax[1]=max(y);

		if(is_test_data==1){
			List ret(13);
			ret[0] = overall_lik;
			ret[1] = overall_trees;
			ret[2] = overall_mat;
			ret[3] = prev_round_preds;
			ret[4] = prev_round_BIC;
			ret[5] = orig_preds;
			ret[6] = overall_overall_sum_trees;
			ret[7] = overall_overall_sum_trees_mat;
			ret[8] = orig_test_preds;
			ret[9] = minmax;
			ret[10] = overall_overall_sum_BIC;
			ret[11] = end_BIC;
			ret[12] = overall_overall_sum_tree_resids;

			return(ret);
		}else{
			List ret(12);
			ret[0] = overall_lik;
			ret[1] = overall_trees;
			ret[2] = overall_mat;
			ret[3] = prev_round_preds;
			ret[4] = prev_round_BIC;
			ret[5] = orig_preds;
			ret[6] = overall_overall_sum_trees;
			ret[7] = overall_overall_sum_trees_mat;
			ret[8] = minmax;
			ret[9] = overall_overall_sum_BIC;
			ret[10] = end_BIC;
			ret[11] = overall_overall_sum_tree_resids;
			return(ret);
		}
	}
}
