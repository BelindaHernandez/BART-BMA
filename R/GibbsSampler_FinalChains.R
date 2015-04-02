#this function will create a final MCMC chain where draws from each sum of tree model in OW are randomly chosen from their
#after burnin iterations. The number of draws included from each sum of tree model is proportional to its posterior probability.
get_final_chain<-function(bartbma,gibbs_output,num_iter,burnin){
  
sum_of_tree_BIC<-unlist(bartbma)
weights<-exp(sum_of_tree_BIC-(max(sum_of_tree_BIC)+log(sum(exp(sum_of_tree_BIC-max(sum_of_tree_BIC))))))
sigma_chains<-gibbs_output[[5]]
final_length=num_iter-burnin
num_its_to_sample<-round(weights*final_length)
final_sigma_chain<-numeric(0)
y_posterior_sum_trees<-gibbs_output[[4]]
final_y_chain<-matrix(nrow=0,ncol=ncol(y_posterior_sum_trees[[1]]))
for(i in 1:length(sigma_chains)){
  sample_its<-sample(burnin:num_iter,num_its_to_sample[i])
  final_sigma_chain<-c(final_sigma_chain,sigma_chains[[i]][sample_its])
  #now do the same for predicted response updates
  post_y_i<-y_posterior_sum_trees[[i]]
  final_y_chain<-rbind(final_y_chain,post_y_i[sample_its,])
}
ret<-list()
length(ret)<-2
ret[[1]]<- final_y_chain
ret[[2]]<-final_sigma_chain
return(ret)
}
