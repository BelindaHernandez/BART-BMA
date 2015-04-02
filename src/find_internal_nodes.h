#include <Rcpp.h>
using namespace Rcpp;
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
