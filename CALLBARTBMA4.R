BART_BMA4=function(x.train,y.train,a=3,nu=3,sigquant=0.9,c=1000,
                   pen=12,num_cp=20,x.test=matrix(0.0,0,0),num_rounds=5,alpha=0.95,beta=1,split_rule_node=0,gridpoint=0,maxOWsize=500
){
  binary=FALSE
  start_mean=0
  start_sd=1
  mu=0
  sigma_mu=0; 
  sigma=sd(y.train)/(max(y.train)-min(y.train))
  qchi = qchisq(1.0-sigquant,nu,1,0);
  lambda = (sigma*sigma*qchi)/nu;
  # print(c("lambda is ",lambda))
  #  lambda=0.1
  if(is.factor(y.train)) {
    # if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    binary = TRUE
    #  y.train = as.numeric(y.train)-1
    stop("Response must be a numeric vector")
  } else {
    if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
      # cat('NOTE: assumming numeric response is binary\n')
      #  binary = TRUE
      stop("Response must be a numeric vector")
    }
  }
  
  if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.train)) x.test = as.matrix(x.test)
  
  if(is.matrix(x.train)) {
    if(nrow(x.test)) {
      if(!is.matrix(x.test)) stop('x.test must be a matrix')
    } 
  }
  
  
  #check input arguments:
  # if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
  #if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
  if((!is.matrix(x.train))) stop("argument x.train must be a double matrix")
  if((!is.matrix(x.test)) ) stop("argument x.test must be a double matrix")
  if(!binary) {
    if((!is.vector(y.train))) stop("argument y.train must be a double vector")
  }
  if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
  if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
  #if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
  #if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
  #if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
  if(c<1)stop("Value of Occam's Window has to be greater than 0."); 
  if(num_cp<0 || num_cp>100)stop("Value of num_cp should be a value between 1 and 100."); 
  
  # sourceCpp("010914_BART_BMA_test_data_PELTmeanvar_scaledresponse6.cpp")  
  start=proc.time()
  test_BMA34=BART_BMA(x.train,y.train,start_mean,start_sd,a,mu,nu,lambda,c,sigma_mu,
                      pen,num_cp,x.test,num_rounds,alpha,beta,split_rule_node,gridpoint,maxOWsize)
  end=proc.time()
  print(c("Time is " ,start-end))
  return(test_BMA34)
}
