library(MASS)
library(msm)
library(coda)
library(pROC)
library(deSolve)
library(gplots)
library(pracma)
library(PerfMeas)
library(texmex)
Rcpp::sourceCpp('~/llFunc.cpp')
#######################################################################################################################
# functions of likelihood,prior of W and parameter estimation
#
# density of observations -----------------------------

mu_1 = 0.5
sigma_1 = 0.01

rate_0=10

ddens_1 = function(j){dlnorm(x = j, meanlog = mu_1, sdlog = sigma_1)}
ddens_0 = function(j){dexp(x = j, rate = rate_0)}



rdens_1 = function(){rlnorm(n = 1, meanlog = mu_1, sdlog =  sigma_1)}
rdens_0 = function(){rexp(n = 1, rate = rate_0)}

#======================================================

# compute dynamical system state-------------------
# W = non-negative and negative adjacency matrix
# q = perturbation vector (1 = perturbed, 0 = unperturbed)
# r0 = initial state of unperturbed dynamical system
# times = times at where to compute system state
# method = which method to use
# returns vector of states for each time point (i.e. a #states x #time points matrix)
compute_states = function(W, q, r0, times, method="ODE",P_S,tstp=tstep){
  method = match.arg(method, several.ok=FALSE)
  if(method == "ODE"){ # linear ODE approach ==> matrix exponential
    nstates=length(r0)
    
    rq=matrix(0,nrow = times,ncol = nstates,dimnames = list(paste0('T_',as.character(1:times)),paste0(colnames(W)[1:nstates])))
    # for(myt in 1:(times-1))
    #  rq[myt+1,] =expAtv(t(W), (c(rq[1,],q)), t=myt)$eAtv[1:nstates]
    rq=ode(y = r0,times = seq(1,times,tstp),func = rfunc, parms = list(W,q,0,P_S))[,-1]  # parms[[3]]=0 in case of no depletion  
    
  }
  else if(method == "difference"){ # linear difference equation
    times = 1:length(times) # just indices of time steps
    
  }
  else if(method == "steady state"){
    
  }
  return(rq)
}
#=======================================================


rfunc <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    
    # import <- input(t)
    #diag(parms[[1]])= diag(parms[[1]])*parms[[3]]
    #dx= t(parms[[1]]) %*% (x) + (c(import,import,import) * parms[[2]])
    dx= (t(parms[[1]]) %*% (x) + (t(parms[[4]])%*% parms[[2]]))
    
    list(c(dx))
  })
}

#=======================================================

# compute log-prior----------------------------
# W = wights of hidden structure
# rho, lambda, nu, tau = hyperparameters of prior
# returns log-prior

prior = function(W, rho, lambda, nu, tau){
  nstates=dim(rho)[1]
  logp = sum(log(( rho*dlnorm(abs(W)[1:nstates,1:nstates],meanlog = nu,sdlog = tau))+ ((1-rho)*dexp(abs(W)[1:nstates,1:nstates],rate = lambda))))
  
  return(logp)
}

#===========================================================

#mixture coefficient of measurement distribution --------------------------------
gammas = function(R,r0#Parameter of network
){
  return(2* (exp(abs(R-r0))/(1+exp(abs(R-r0))))-1)
  
}

#========================================================

# learn a model from data-------------------
# dat = see above
# Q = matrix of perturbation vectors (each column = one perturbation vector)
# rho = optional matrix of prior edge probabilities or just a constant
# hyper.W = optional list of additional hyper-parameters for the W-prior (lambda, nu, tau)
# prior.theta = optional matrix of prior attachement probabilities of E-nodes to S-nodes
# seed = random seed
# returns model structure
learn = function(dat, Q, rho, hyper.W, prior.theta,seed=123,P_S,ts){
  npert = NROW(Q)
  nstates = NCOL(Q)
  if(is.null(rho))
    rho = matrix(0.5, ncol=nstates, nrow=nstates,	dimnames= list( colnames(W), colnames(W)))
  if(is.null(prior.theta))
    prior.theta = matrix(1/(nstates+1), nrow=NROW(dat), ncol=nstates)
  if(is.null(hyper.W)){
    hyper.W$lambda = 10 # scale parameter of exponential distribution
    hyper.W$nu = hyper.W$tau = 10 # shape and rate parameters of gamma distribution
  }
  
  times =dim(dat)[3] #tryCatch(as.numeric(dimnames(dat)[[3]]), error = function(e) stop("dimnames[[3]] has to indicate measurement time points"))
  # set.seed(seed)
  
  #W=null_model$W#######################################################################
  
  model = list(W=W, Q=Q, prior.theta=prior.theta, prior.W=list(rho=rho, lamdba=hyper.W$lambda, nu=hyper.W$nu, tau=hyper.W$tau)
               , orig.dat=dat, r0=rep(0, nstates),PS=P_S) # orig.dat should be changed 
  
  
  
  
  
  
  pre_score = prior(model$W,model$prior.W$rho,model$prior.W$lamdba,model$prior.W$nu,model$prior.W$tau)+Likelihood(model,dat,compute_states ,gammas ,ddens_1,ddens_0,P_S,ts)
  
  # perform MCMC
  acc_sd=list()
  rej_sd=list()  
  acc_res=list()
  acc_ratio=list()  
  rej_ratio=list()
  rej_res=list()
  all_res=list()
  all_ratio=list()
  all_flag=list()
  all_Q=list()
  #suggested_sds = 2^((-10):1)
  #suggested_sds = seq(2^(-10),0.2,length=5000)
  for(itr in 1:100)
  { 
    
  Q=0
  move_type = sample(x = c(1,2,3,4,5),size = 1,prob = c(1/5,1/5,1/5,1/5,1/5))
  
  indx = sample(1:nstates,2,replace = TRUE)
  i_indx = indx[1]
  j_indx= indx[2]
  
  
  new_model=model
  
  flag_move=0
  
  if(move_type == 1)      #insert new edge
  {
    if((model$W[i_indx,j_indx]==0)){ # no bidirectional edge !!! let loops
      
      
      flag_move = 1
      sign_diag <-1
      if(i_indx== j_indx)
        sign_diag <- sample(c(1,-1),1)
      
      new_model$W[i_indx,j_indx] = sign_diag * exp( rnorm(1, mean = new_model$prior.W$nu, sd = new_model$prior.W$tau)) # new_w_ij = lnorm(nu, tau)
      
      Q_new_cond_old <- dlnorm(abs(new_model$W[i_indx,j_indx]), meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/5) * (1/sum(model$W==0))#insert to old
      
      Q_old_cond_new <- (1/5) * (1/sum(new_model$W!=0)) #delete from new
      
      if(model$W[j_indx,i_indx]!=0 & i_indx!=j_indx){
        Q_new_cond_old <- Q_new_cond_old + (dlnorm(new_model$W[i_indx,j_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/5) * (1/sum(xor(c(model$W[lower.tri(model$W)]),c(model$W[upper.tri(model$W)]))))) #reverse from old
        
      }
      
      if(new_model$W[j_indx,i_indx]==0 & i_indx!=j_indx){
        Q_old_cond_new <- Q_old_cond_new + ( (1/5) * (1/sum(xor(c(new_model$W[lower.tri(new_model$W)]),c(new_model$W[upper.tri(new_model$W)]))))) #reverse from new
        
      }
      if(Q_new_cond_old&&Q_old_cond_new)
        
        Q<- log(Q_old_cond_new)- log(Q_new_cond_old)
    }
    
    # if (!is.DAG(new_model$W)){#reject new_model if its DAG
    #   flag_move = 0
    #   Q = 0
    # }
  } else if(move_type==2)#delete an edge
  {
    if(model$W[i_indx,j_indx]!=0){
      
      flag_move = 1
      new_model$W[i_indx,j_indx] = 0
      
      Q_new_cond_old <-  (1/5) * (1/sum(model$W!=0))#DELET FROM OLD
      
      Q_old_cond_new <- dlnorm(abs(model$W[i_indx,j_indx]), meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) *(1/5) * (1/sum(new_model$W==0)) #insert to new
      
      if(model$W[j_indx,i_indx]==0 & i_indx!=j_indx){
        Q_new_cond_old <- Q_new_cond_old + ( (1/5) * (1/sum(xor(c(model$W[lower.tri(model$W)]),c(model$W[upper.tri(model$W)]))))) #reverse from old
        
      }
      
      if(new_model$W[j_indx,i_indx]!=0 & i_indx!=j_indx){
        Q_old_cond_new <- Q_old_cond_new + ( dlnorm(model$W[i_indx,j_indx], meanlog = model$prior.W$nu, sdlog = model$prior.W$tau) *(1/5) * (1/sum(xor(c(new_model$W[lower.tri(new_model$W)]),c(new_model$W[upper.tri(new_model$W)]))))) #reverse from new
        
      }
      if(Q_new_cond_old&&Q_old_cond_new)
        
        Q<- log(Q_old_cond_new)- log(Q_new_cond_old)
    }
    
  }else if(move_type==3)#reverse an edge
  {
    if(i_indx!=j_indx){
      if( xor( model$W[i_indx,j_indx], model$W[j_indx,i_indx]) ){
        
        flag_move = 1
        new_model$W[i_indx,j_indx]=model$W[j_indx,i_indx]
        new_model$W[j_indx,i_indx]=model$W[i_indx,j_indx]
        
        if(new_model$W[i_indx,j_indx]!= 0){
          
          Q_new_cond_old <- (1/5) * dlnorm(new_model$W[i_indx,j_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/sum(xor(c(model$W[lower.tri(model$W)]),c(model$W[upper.tri(model$W)]))))   #reverse from old
          
          Q_old_cond_new <- (1/5) * dlnorm(new_model$W[i_indx,j_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/sum(xor(c(new_model$W[lower.tri(new_model$W)]),c(new_model$W[upper.tri(new_model$W)]))))#reverse from new
          
        } else {
          
          Q_new_cond_old <- (1/5) * dlnorm(new_model$W[j_indx,i_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/sum(xor(c(model$W[lower.tri(model$W)]),c(model$W[upper.tri(model$W)]))))#reverse from old
          
          Q_old_cond_new <- (1/5) * dlnorm(new_model$W[j_indx,i_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) * (1/sum(xor(c(new_model$W[lower.tri(new_model$W)]),c(new_model$W[upper.tri(new_model$W)]))))#reverse from new
          
        }
        if(Q_new_cond_old&&Q_old_cond_new)
          
          Q<- log(Q_old_cond_new)- log(Q_new_cond_old)  # should be one
        
        
        
      }}
  } else if(move_type==4)#swip weights
  {
    if(i_indx!=j_indx){
      if(model$W[i_indx,j_indx]!= 0) {
        
        pair<-expand.grid(1:nstates,1:nstates)
        pair <- pair[-which(pair[,1]==pair[,2]),]
        indx <-sample(1:nrow(pair),1)
        k_indx <- pair[indx,1]
        l_indx <- pair[indx,2]
        
        if((i_indx!=k_indx)|(j_indx!=l_indx)){
          if((model$W[i_indx,j_indx]!= model$W[k_indx,l_indx]) &  (model$W[k_indx,l_indx] != 0))
          {
            flag_move = 1
            
            new_model$W[k_indx,l_indx]=model$W[i_indx,j_indx]
            new_model$W[i_indx,j_indx]=model$W[k_indx,l_indx]
            
            Q_new_cond_old =  1/5 #(1/sum(c(new_model$W[lower.tri(new_model$W)|upper.tri(new_model$W)])==model$W[i_indx,j_indx]))  dlnorm(new_model$W[i_indx,j_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) *(1/5) * (1/sum(c(new_model$W[lower.tri(new_model$W)|upper.tri(new_model$W)])!=0)) * unique
            
            Q_old_cond_new = 1/5 #dlnorm(model$W[i_indx,j_indx], meanlog = model$prior.W$nu, sdlog = model$prior.W$tau) *(1/5) * (1/sum(c(new_model$W[lower.tri(new_model$W)|upper.tri(new_model$W)])!=0)) * unique
            
            Q<- log(Q_old_cond_new)- log(Q_new_cond_old)
            
          }
        }else{
          if((model$W[i_indx,j_indx]!= model$W[l_indx,k_indx]) &  (model$W[l_indx,k_indx] != 0))
          {
            flag_move = 1
            
            new_model$W[l_indx,k_indx]=model$W[i_indx,j_indx]
            new_model$W[i_indx,j_indx]=model$W[l_indx,k_indx]
            
            Q_new_cond_old =   1/5 #dlnorm(new_model$W[i_indx,j_indx], meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) *(1/5) * (1/sum(c(new_model$W[lower.tri(new_model$W)|upper.tri(new_model$W)])!=0)) * unique
            
            Q_old_cond_new =  1/5 #dlnorm(model$W[i_indx,j_indx], meanlog = model$prior.W$nu, sdlog = model$prior.W$tau) *(1/5) * (1/sum(c(new_model$W[lower.tri(new_model$W)|upper.tri(new_model$W)])!=0)) * unique
            
            Q<- log(Q_old_cond_new)- log(Q_new_cond_old)
          }
        }
        
      }}
    
  } else #add value to edge
  {
    if(model$W[i_indx,j_indx]!=0){
      
      flag_move = 1
      sign_diag <-1
      if(i_indx== j_indx)
        sign_diag <- sample(c(1,-1),1)
      new_model$W[i_indx,j_indx] = sign_diag * exp(log(abs(model$W[i_indx,j_indx]))+rnorm(1,mean = 0,sd = 0.02))
      if (new_model$W[i_indx,j_indx]!= 0){
        
        Q_new_cond_old = dlnorm(abs(new_model$W[i_indx,j_indx]), meanlog = new_model$prior.W$nu, sdlog = new_model$prior.W$tau) #*(1/5) * (1/sum(model$W!=0))
        
        Q_old_cond_new = dlnorm(abs(model$W[i_indx,j_indx]), meanlog = model$prior.W$nu, sdlog = model$prior.W$tau) #* (1/5) * (1/sum(model$W!=0))
        
        if(Q_new_cond_old&&Q_old_cond_new)
          Q<- log(Q_old_cond_new)- log(Q_new_cond_old)
        
      }else{
        
        flag_move = 0
        Q = 0
        
      }
    }
    
  }
  
  
  l=Likelihood(new_model,dat,compute_states ,gammas ,ddens_1,ddens_0,PS,ts)
  if(Q){
    new_score = prior(new_model$W,new_model$prior.W$rho,new_model$prior.W$lamdba,new_model$prior.W$nu,new_model$prior.W$tau)+l+log(Q_old_cond_new)
    pre_score = prior(model$W,model$prior.W$rho,model$prior.W$lamdba,model$prior.W$nu,model$prior.W$tau)+Likelihood(model,dat,compute_states ,gammas ,ddens_1,ddens_0,PS,ts)+ log(Q_new_cond_old)
  }else{
    
    new_score = prior(new_model$W,new_model$prior.W$rho,new_model$prior.W$lamdba,new_model$prior.W$nu,new_model$prior.W$tau)+l
    pre_score = prior(model$W,model$prior.W$rho,model$prior.W$lamdba,model$prior.W$nu,model$prior.W$tau)+Likelihood(model,dat,compute_states ,gammas ,ddens_1,ddens_0,PS,ts) 
  }
  all_res=c(all_res,list(new_model$W))
  all_ratio=c(all_ratio,l)
  
  # logratio= new_score - pre_score + Q 
  randratio = (runif(1,min=0,max=2))
  if (((new_score>pre_score)&&(flag_move==1))||((exp(new_score-pre_score)>=randratio) &&(flag_move==1)))
  {
    print("itr")
    print(Q)
    print(move_type)
    
    model=new_model
    acc_res= c(acc_res,list(new_model$W))
    acc_ratio=c(acc_ratio,list(l))
    all_flag=c(all_flag,1)
    pre_score = new_score
    all_Q = c(all_Q,c(Q,move_type))
  } else {
    rej_res= c(rej_res,list(new_model$W))
    rej_ratio=c(rej_ratio,list(l))
    all_flag=c(all_flag,0)
    all_Q = c(all_Q, c(Q,move_type))
  }
  
}

  return(list(all_res,all_ratio,all_flag,acc_res,acc_ratio,rej_res,rej_ratio,all_Q))
}
#========================================================

# compute log-likelihood according to formula in supplements---------------------------------
# model = list(W, Q, prior.theta, prior.W=list(rho, lambda, nu, tau), orig.dat, r0)
# dat = data array (#E-nodes x #perturbations x time points)
# returns log-likelihood
likelihood = function(model, dat){
  
  nstates = NCOL(model$W)   #S numbers
  Snames <- colnames(model$W)
  Pnames <- colnames(model$Q)
  Enames <- rownames(dat)[-dim(dat)[[1]]]
  times = as.numeric(dim(dat)[[3]]) 
  ll = 0
  R=list()
  gamma_val=list()
  
  for(q in 1:NROW(model$Q)){ 
    
    # iterate over all perturbation experiments
    R[[q]] = compute_states(model$W, model$Q[q,], model$r0, times,"ODE",model$PS)
    temp=R[[q]][which(abs(R[[q]])>700,arr.ind = T)]                #thereshold over states because of exp() in gamma function
    R[[q]][which(abs(R[[q]])>700,arr.ind = T)]=700*sign(temp)
    
    
    # compute gammas
    gamma_val[[q]]=gammas(R[[q]],model$r0 )
    # compute Eq. (3), ensuring you never take log(0)
    # plug into formula in supplements:
    for(i in 1:(dim(dat)[1]-1)){ 
      # iterate over all E-nodes
      # add the second formula in the Supplements to ll
      temp=c()
      for(s in 1:nstates)
      {
        temp[s]<- sum(log((((gamma_val[[q]][,s]*ddens_1(dat[i,q,],q)) + ((1-gamma_val[[q]][,s])*ddens_0(dat[i,q,],q)))*(model$prior.theta[i,s]))+1))
      }
      
      a = max(temp)
      small_a=log(sum(exp(temp-a)))
      #track likelihood computation for each E_node and different parents
      #print(a)
      #print(temp)
      #print(temp-a)
      #print(small_a)
      ll=ll+sum(a,small_a)
    }
    
  }
  return(ll)  
}
#=======================================================

# analyse and plot score of MCMC results---------------
plotMCMC<-function(model, results){
  

  plot(x,y,type = "l")
  
}

#############################################################
postProcess<- function(results,model){
  
  nstates = ncol(model$W)
  
  len=length(results)
  b_in=len/2
  thn=10
  
  postPrior=matrix(0,nstates,nstates)
  sc=matrix(0,nstates,nstates)
  med=matrix(0,nstates,nstates)
  
  res=list()
  
  for(from in 1:nstates)
    for(to in 1:nstates)
    {
      if(from!=to)
      {
        W_list=sapply(1:len,function(i) unlist(results[[i]])[from,to])
        
        postPrior[from,to]=sum(Thin(W_list[b_in:len],By = thn)!=0)/(length(Thin(W_list[b_in:len],By = thn)))
        
        if(sd(Thin(W_list[b_in:len],By = thn))){
          sc[from,to] =mean(Thin(W_list[b_in:len],By = thn))/sd(Thin(W_list[b_in:len],By = thn))
        } else{
          sc[from,to] =   mean(Thin(W_list[b_in:len],By = thn))
        }
        
        med[from,to] = median(Thin(W_list[b_in:len],By = thn))/median(abs(Thin(W_list[b_in:len],By = thn)-median(Thin(W_list[b_in:len],By = thn))))
        
      }
    }
  
  return(list(postPrior,sc))}
#############################################################
ROC_curve<- function(postPrior, sc,model){
  
  nstates = ncol(model$Q)
  model$W= model$W[1:nstates,1:nstates]
  
  #pROC
  x=roc( response =   as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)])),predictor =     as.vector(c(postPrior[lower.tri(postPrior)],postPrior[upper.tri(postPrior)])),auc = TRUE,plot = TRUE)
  print(auc(x))
  #ifile <- paste0(sim_num,'pprior_pROC_ROC.pdf')
  #pdf(ifile)
  #plot.roc(x,legacy.axes = TRUE,print.auc = TRUE,auc.polygon = TRUE)
  #d <- dev.off()
  
  y=roc( response =   as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)])),predictor =   as.vector(c(sc[lower.tri(sc)],sc[upper.tri(sc)])),auc = TRUE,plot = TRUE)
  print(auc(y))
  #ifile <- paste0(sim_num,'sc_pROC_ROC.pdf')
  #pdf(ifile)
  #plot.roc(y,legacy.axes = TRUE,print.auc = TRUE,auc.polygon = TRUE)
  #d <- dev.off()
  write.table(c(x$auc), file = paste0('Res_ROC_posteriorP'),append = TRUE, sep = ",",row.names = FALSE,col.names = FALSE)
  write.table(c(y$auc), file =  paste0('Res_ROC_Score'),append = TRUE, sep = ",",row.names = FALSE,col.names = FALSE)
  
}
#############################################################
PR_curve<- function(postPrior, sc ,model){
  
  nstates = ncol(model$Q)
  model$W= model$W[1:nstates,1:nstates]
  
  #prefMeas
  #print(AUC.single(pred = as.vector(c(postPrior[lower.tri(postPrior)],postPrior[upper.tri(postPrior)])),labels = as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)]))))
  #print(AUC.single(pred = as.vector(c(sc[lower.tri(sc)],sc[upper.tri(sc)])),labels = as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)]))))
  res=list(precision.at.all.recall.levels(scores = as.vector(c(postPrior[lower.tri(postPrior)],postPrior[upper.tri(postPrior)])), labels = as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)]))))
  AUPRC_p <- AUPRC (res, comp.precision=TRUE)
  print(AUPRC_p)
  #ifile <- paste0(sim_num,'pprior_prefMeas_pRec.pdf')
  #pdf(ifile)
  #precision.recall.curves.plot(res,plot.precision = TRUE)
  #d <- dev.off()
  precision.recall.curves.plot(res,plot.precision = TRUE)
  
  res=list(precision.at.all.recall.levels(scores = as.vector(c(sc[lower.tri(sc)],sc[upper.tri(sc)])), labels = as.vector(c(model$W[lower.tri(model$W)],model$W[upper.tri(model$W)]))))
  AUPRC_s <- AUPRC (res, comp.precision=TRUE)
  print(AUPRC_s)
  #ifile <- paste0(sim_num,'sc_prefMeas_pRec.pdf')
  #pdf(ifile)
  #precision.recall.curves.plot(res,plot.precision = TRUE)
  #d <- dev.off()
  precision.recall.curves.plot(res,plot.precision = TRUE)
  
  write.table(c(AUPRC_p), file =  paste0('Res_PR_posteriorP'),append = TRUE, sep = ",",row.names = FALSE,col.names = FALSE)
  write.table(c(AUPRC_s), file =  paste0('Res_PR_Score'),append = TRUE, sep = ",",row.names = FALSE,col.names = FALSE)
  
}
#############################################################
convergence = function(results,snum,itr){
  
  len=min(sapply(1:itr,function(i) length(results[[i]][[4]])))
  b_in=len/2
  
  if(itr>1){
    for(i in 1:(itr-1)){
      combinedchains = mcmc.list(mcmc(data = unlist(results[[i]][[5]][b_in:len])),mcmc(unlist(results[[i+1]][[5]] [b_in:len])))
      ifile <- paste0(i,i+1,'combinedchains.pdf')
      pdf(ifile)
      par(mfrow=c(2,2),mar=c(3,4,1,1))
      plot(combinedchains,main=paste0(i,i+1,'all_scoreCombined'))
      g=gelman.diag(combinedchains)
      gelman.plot(combinedchains,main=paste0('scoreCombined'))
      save(g,file = paste0('gelman_all.RData'))
      d <- dev.off() 
      
      for(from in 1:snum)
      {
        ifile <- paste0(i,i+1,from,'combinedchains.pdf')
        pdf(ifile)
        par(mfrow=c(2,2),mar=c(3,4,1,1))
        
        for(to in 1:snum){
          if(to != from){
            combinedchains = mcmc.list(mcmc(sapply(b_in:len,function(j) unlist(results[[i]][[4]][[j]])[from,to])),mcmc(sapply(b_in:len,function(j) unlist(results[[i+1]][[4]][[j]])[from,to])))
            plot(combinedchains,main=paste0(from,'---->',to,'W_ijCombined'))
            g=gelman.diag(combinedchains)
            gelman.plot(combinedchains,main=paste0(from,'---->',to,'W_ijCombined'))
            print(gelman.diag(combinedchains))
            save(g,file = paste0(i,'compared',i+1,'from',from,'to',to,'gelman.RData'))
          }}
        d <- dev.off() 
      }
    }}
}

#########################################################
