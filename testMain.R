source("~/testMCMC.R")


S_num<-3
E_num<-10
times<-10
tstep <-0.1

#make test true W
W<-matrix(0,3,3)
W[1,3]<-1
W[2,3]<-1
diag(W)=-0.1

PS<- diag(1,3) #P to S connections


l <- rep(list(0:1), S_num)
Q = as.matrix(expand.grid(l))    # known state of intervention nodes (interventing one Pnode at each experiment)

O = array(0, dim=c(E_num*S_num, nrow(Q), length(seq(1,times,tstep))), dimnames=list(paste0('E_',as.character(1:(E_num*S_num))), paste0('pert_',as.character(1:NROW(Q))), paste0('T_',1:(length(seq(1,times,tstep)))))) # prepare observation data

X=diag(x = 1,nrow =S_num )
prior.theta = do.call(rbind,lapply(1:S_num ,function(s)rbind(matrix(X[s,],nrow=E_num,ncol=S_num ,byrow=TRUE))))  #S_E connections


r0<-c(0,0,0)  #start state


for(q in 1:nrow(Q)){ 
  
  #propagation under each perturbation
  rq= compute_states(W = W,q = Q[q,],r0 = r0,times = times,method = "ODE",P_S = PS,tstp = tstep)
  
  prob <- gammas(R = rq,r0 = r0)
  
  
  plot(prob[,1],type = "l")
  points(prob[,2],col="red",type = "l")
  points(prob[,3],col="blue",type="l")
  
  
  
  plot(rq[,1])
  points(rq[,2],col="red")
  points(rq[,3],col="blue")
  
  
  
  #simulate observation
  for(i in 1:dim(O)[1]){ # iterate over all E-nodes
    for(t in 1:(length(seq(1,times,tstep)))){# iterate over time points
      
      # draw O[i,q,t] from mixture distribution (Eq. 3)
      indicator = sample(c(FALSE,TRUE),size = 1,prob=c(1-prob[t,which(prior.theta[i,]==1)],prob[t,which(prior.theta[i,]==1)]))
      O[i,q,t]= ifelse(indicator,rdens_1(),abs(rdens_0()))

    }
  }
}


#observation heatmaps
my_palette <- colorRampPalette(c("blue","gray","red"))
for (i in 1:nrow(Q)) {
  heatmap.2((O[,i,]),density.info = "none",col=my_palette,trace="none",Rowv = "none",Colv = "none",scale="none",main=paste("Exp_",i))
  
}
#d <- dev.off()


#truth model                                   
true_model = list(W=W, Q=Q, prior.theta=prior.theta, prior.W=NULL, orig.dat=O, r0=r0, PS=PS)
true_model$prior.W=list(rho=NULL, lambda=NULL, nu=NULL, tau=NULL)  #weight matrix
true_model$prior.W$lambda=1   #parameter for exponential function(nonedge)
true_model$prior.W$nu=-0.2        #parameter for lnorm function(edge)lmean
true_model$prior.W$tau=0.5      #parameter for lnorm function(edge)lsd

true_ll<-Likelihood(true_model,true_model$orig.dat,compute_states ,gammas,ddens_1,ddens_0,true_model$PS,tstep)
##################run mcmc and plots##############################################

res_mcmc <-learn(O, Q, rho=NULL, hyper.W=true_model$prior.W, prior.theta =true_model$prior.theta,seed=123,P_S = true_model$PS,tstep)


plot(unlist(res_mcmc[[2]])) #all_ll
points(unlist(res_mcmc[[5]]),col="red") #accepted_ll
points(unlist(res_mcmc[[7]]),col="blue") #rejected_ll
points(true_ll,col="green") #true_ll


plot(x=c(1:length(res_mcmc[[5]]), y=unlist(res_mcmc)),type = "l") #trace of accepted_ll



plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][1,1]}),type = "l") #trace of acc w for each edge
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][1,2]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][1,3]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][2,1]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][2,2]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][2,3]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][3,1]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][3,2]}),type = "l")
plot(1:length(res_mcmc[[5]]),sapply(1:length(res_mcmc[[5]]), function(i){res_mcmc[[4]][[i]][3,3]}),type = "l")
                                   
                                   
                                   

postMatrices = postProcess(res_mcmc[[4]],true_model)


ROC_curve(postPrior = postMatrices[[1]],sc = postMatrices[[2]] ,model = true_model)


PR_curve(postMatrices[[1]],postMatrices[[2]] ,true_model)  

