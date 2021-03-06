source("~/test/testMCMC.R")
library(igraph)

#######################################################################################################################
# functions for sampling network and data, and modeling perturbation
#######################################################################################################################

#######################################################################################################################
# functions for sampling network and data, and modeling perturbation
#######################################################################################################################

#construct a network structure ----------------------------------
samplenetwork<- function (S_num,P_num,E_num,theta,stype){
  
  #Degree sequences for network motifs
  if(grepl(stype,"FF")){
    deg_seq=c(2,1,0) # You can fix (out)degree of nodes here 
    toedge = paste("S",c(2,3,3),sep = "")
  } else if(grepl(stype,"FB")){
    deg_seq=c(1,1,1) # You can fix (out)degree of nodes here
    toedge = paste("S",c(2,3,1),sep = "")
  } else if(grepl(stype,"Dia")){
    deg_seq=c(2,1,1,0) # You can fix (out)degree of nodes here 
    toedge = paste("S",c(2,3,4,4),sep = "")
  } else if(grepl(stype,"Bi")){
    deg_seq=c(2,0,0,2) # You can fix (out)degree of nodes here 
    toedge = paste("S",c(2,3,2,3),sep = "")
  }  else if(grepl(stype,"PC")){
    deg_seq=c(1,1,1,1)} # You can fix (out)degree of nodes here 
  # else if(grepl(stype,"JK")){
  #   deg_seq=c(1,1,1,0) # You can fix (out)degree of nodes here 
  # }else if(grepl(stype,"UVB")){
  #   deg_seq=c(4,4,3,3,4,1,2,1,1,2,3) # You can fix (out)degree of nodes here 
  # }
  
  newgraph= graph(NULL,n=S_num+P_num+(E_num*S_num),directed = T)    # make an empty graph from S_nodes, E_nodes and P_nodes
  S_names=paste("S",1:S_num,sep="")                         #set nodes' names
  P_names=paste("P",1:P_num,sep="")
  E_names=paste("E",1:(E_num*S_num),sep="")
  rownames(theta)=E_names
  V(newgraph)$name<-c(S_names,P_names,E_names)              
  V(newgraph)$color<-c( paste(rep("yellow",S_num),sep=""), paste(rep("red",P_num),sep=""), paste(rep("white",(E_num*S_num)),sep=""))             #set node's color
  
  
  fromedge = rep(S_names,deg_seq)         
  # toedge = NULL
  # for (j in 1:S_num){
  #   if(stype %in% c("FF","FB","PC")){
  #     toedge = c(toedge, sample(S_names[1:S_num],deg_seq[j]))               # for FF,FB,PC structure
  #   }else{toedge = c(toedge, sample(S_names[(j+1):S_num],deg_seq[j]))}           #sample edges between S_S nodes   
  #   if(grepl(stype,"UVB")){
  #     toedge = NULL
  #     toedge = paste("S",c(2,3,4,5,1,3,10,11,1,2,9,5,6,7,1,4,7,9,4,4,8,6,1,2,11,2,9,10),sep = "")
  #   }
  # }
  # 
  
  edgelist = rbind(fromedge,toedge)
  newgraph = add_edges(newgraph,edgelist)        #add S_S edges to the empty graph
  
  
  #Sample scale free network for S nodes---------------------------
  ##! degs <- sort(sample(1:S_num, S_num, replace=TRUE, prob=(1:S_num)^-2),decreasing=TRUE)     #generate Power-law degree distribution for Snodes graph
  ##!rand_graph <- sample_degseq(degs,rev(degs),method = "simple.no.multiple",)         #use pawer-law degree for generating scale free network  
  
  
  
  #connect P layer and S layer 
  edgelist = (rbind(P_names,S_names))      
  newgraph = add_edges(newgraph,edgelist)        #add P_S edges to the graph          
  
  
  #connect S layer and E layer  
  # if (E_num < S_num) stop("Number of E-nodes is smaller than the number of signal nodes. The model is unidentifiable.")
  edgelist =  as.vector(rbind(sapply(1:(E_num*S_num), function(j) S_names[which(theta[E_names[j],]==1)]),E_names))
  newgraph = add_edges(newgraph,edgelist)        #add S_E edges to the graph          
  
  return(newgraph)  
}
#==========================================================

#update hidden nodes state using "ODE" ------------------------------------

uVal= function(t){
  
  return(1-(1/(1+t)))
}



Bifan <- function(t, x,parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1=(U*parms[1])- (parms[5] * x[1])   
    dx2=(U*parms[2])- (parms[5] * x[2])   
    dx3=(U*parms[3])+x[2]+x[1]- (parms[5] * x[3])   
    dx4=(U*parms[4])+x[2]+x[1]- (parms[5] * x[4])   
    list(c(dx1,dx2,dx3,dx4))
  })}

Diamond <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1=(U*parms[1])- (parms[5] * x[1])   
    dx2=(U*parms[2])+x[1]- (parms[5] * x[2])   
    dx3=(U*parms[3])+x[1]- (parms[5] * x[3])   
    dx4=(U*parms[4])+x[2]+x[3]- (parms[5] * x[4])   
    list(c(dx1,dx2,dx3,dx4))
  })}

FeedForward <- function(t, x,parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1=(U*parms[1])- (parms[4] * x[1])   
    dx2=(U*parms[2])+x[1]- (parms[4] * x[2])   
    dx3=(U*parms[3])+x[2]+x[1]- (parms[4] * x[3])   
    list(c(dx1,dx2,dx3))
  })}
FeedBackward <- function(t, x,parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1=(U*parms[1])+ x[2]- (parms[4] * x[1])   
    dx2=(U*parms[2])+x[3]- (parms[4] * x[2])   
    dx3=(U*parms[3])+x[1]- (parms[4] * x[3])   
    list(c(dx1,dx2,dx3))
  })}


ProteinCascade <- function(t, x,parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1=(U*parms[1])- (parms[5] * x[1])   
    dx2=(U*parms[2])+(x[1]/(1+x[4]))- (parms[5] * x[2])   
    dx3=(U*parms[3])+x[2]- (parms[5] * x[3])   
    dx4=(U*parms[4])+x[3]- (parms[5] * x[4])   
    list(c(dx1,dx2,dx3,dx4))
  })}

JAKSTAT <- function(t, x,parms) {
  with(as.list(c(parms, x)),{
    U=uVal(t)
    dx1= -(parms[1]* x[1]* u)  
    dx2= -(parms[1]* x[1]* u) -2 * (parms[2] * x[2]* x[2])   
    dx3= (parms[2] * x[2]) - (parms[3] * x[3])   
    dx4= (parms[3] * x[3])   
    list(c(dx1,dx2,dx3,dx4))
  })}

UVB <- function(t, x, parms){
  
  dx1 = -2 * ka1 * x[1] * x[1] * x[4] * x[4] + 2 * kd1 * x[5] +
    ks1 * (1 + UV * n3 * (x[9] +  FHY3)) -
    kdr1 * (1 + ( n1 * UV )) * x[1] - kd2 * x[3] -
    2 * ka2 * x[1] * x[1] * x[2]
  
  dx2 = -ka2 * x[1] * x[1] * x[2] + kd2 * x[3] -
    ka4 * x[2] * x[10] + kd4 * x[11]
  
  dx3 = -kd2 * x[3] + ka2 * x[1] * x[1] * x[2]
  
  dx4 = -2 * k1 * x[4] * x[4] + 2 * k2 * UV * x[6] -
    2 * ka1 * x[1] * x[1] * x[4] * x[4] + 2 * kd1 * x[5] -
    ka3 * x[4] * x[7]
  
  dx5 = -kd1 *  x[5] + ka1 * x[1] * x[1] * x[4] * x[4]
  
  dx6 = -k2 * x[6] + k1 * x[4] * x[4] + kd3 * x[8] * x[8]
  
  dx7 = -ka3 * x[4] * x[7] + ks2 * ( 1 + UV * x[5]) -
    kdr2 * x[7] + 2 * kd3 * x[8] * x[8]
  
  dx8 = -2 * kd3 * x[8] * x[8] + ka3 * x[4] * x[7]
  
  dx9 = -kdr3 * (x[3] / (kdr3a + x[3]) + x[11] / (kdr3b + x[11]) ) * x[9] +
    k3sp * ( 1 + n2 * UV ) - kdr3 * (x[5] / (ksr + x[5])) * x[9]
  
  dx10 = -ka4 * x[2] * x[10] + kd4 * x[11]
  
  dx11 = -kd4 * x[11] + ka4 * x[2] * x[10]
}

UVB1 <- function(t, x, parms){
  with(as.list(c(parms, x)),{
    dx1 = 2 *( -ka1 * x[1] * x[1] * x[4] * x[4] + kd1 * x[5] +
                 kd2 * x[3] - ka2 * x[1] * x[1] * x[2])
    
    dx2 = -ka2 * x[1] * x[1] * x[2] + kd2 * x[3] -
      ka4 * x[2] * x[10] + kd4 * x[11]
    
    dx3 = -kd2 * x[3] + ka2 * x[1] * x[1] * x[2]
    
    dx4 = 2 * (-k1 * x[4] * x[4] + k2 * UV * x[6] -
                 ka1 * x[1] * x[1] * x[4] * x[4] + kd1 * x[5])
    
    dx5 = -kd1 *  x[5] + ka1 * x[1] * x[1] * x[4] * x[4]
    
    dx6 = -k2 * UV * x[6] + k1 * x[4] * x[4] 
    
    dx7=0
    dx8=0
    dx9=0
    
    dx10 = -ka4 * x[2] * x[10] + kd4 * x[11]
    dx11 = -kd4 * x[11] + ka4 * x[2] * x[10]
    
    list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10,dx11))
  })}
#===========================================================

#simulate data for a model ---------------------------------------
# model = list(W, Q, prior.theta, prior.W=list(rho, lambda, nu, tau), orig.dat, r0)
# N = number of observations to simulate
# times = times points to simulate
# return data array (N x #states x times)
sim.data = function(model,  nE, times,stype ,ts){
  npert = NROW(model$Q) #number of P nodes
  R=list()
  gamma_val=list()
  O = array(0, dim=c(nE, npert, length(seq(1,times,tstep))), dimnames=list(paste0('E_',as.character(1:nE)), paste0('pert_',as.character(1:npert)), paste0('T_',1:length(seq(1,times,tstep))))) # prepare observation data
 
  for(q in 1:npert){ 
    
    
    
    # draw O[i,q,t] from mixture distribution (Eq. 3)
    if(grepl(stype,"FF")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=FeedForward)[,-1] #no seld-depletion
    } else if(grepl(stype,"FB")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=FeedBackward)[,-1] #no seld-depletion
    }else if(grepl(stype,"Dia")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=Diamond)[,-1] #no seld-depletion
    } else if(grepl(stype,"Bi")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=Bifan)[,-1] #no seld-depletion
    } else if(grepl(stype,"PC")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=ProteinCascade)[,-1] #no seld-depletion
    }else if(grepl(stype,"JK")){
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c(model$Q[q,],0),func=JAKSTAT)[,-1] #no seld-depletion
    }else if(grepl(stype,"UVB")){
      #R[[q]]=ode(y= model$r0, times =seq(0,times,by=1), parms = c( k1 = 0.0043 ,k2 = 161.62 , ka1 = 0.0372, ka2 = 0.0611, ka3 = 0, ka4 = 10.1285, kd1 = 94.3524, kd2 = 50.6973, kd3 = 0 , kd4 = 1.1999, ks1 = 0, ks2 = 0, k3sp = 0, ksr = 0, kdr1 = 0, kdr2 = 0, kdr3 = 0, kdr3a = 0, kdr3b = 0,  UV = 1, FHY3 = 0, n1 = 0, n2 = 0, n3 = 0),func=UVB)[,-1] 
      R[[q]]=ode(y= model$r0, times =seq(0,times,by=ts), parms = c( k1 = 0.0043 ,k2 = 161.62 , ka1 = 0.0372, ka2 = 0.0611, ka4 = 10.1285, kd1 = 94.3524, kd2 = 50.6973 , kd4 = 1.1999, UV = 1),func=UVB1)[,-1]  #UVB1
    }
    
    #dtime= lapply(1:NCOL(model$Q),function(x) lapply(seq(1,(times-1),1), function(z) if ((R[[q]][z+1,x]-R[[q]][z,x])< (0.01*R[[q]][z,x])){z}))
    
    temp=R[[q]][which(abs(R[[q]])>700,arr.ind = T)]            #thereshold over states because of exp() in gamma function
    R[[q]][which(abs(R[[q]])>700,arr.ind = T)]=700*sign(temp)
    
    gamma_val[[q]]=gammas(R[[q]],model$r0 )
    for(i in 1:dim(O)[1]){ # iterate over all E-nodes
      for(t in 1:length(seq(1,times,tstep))){# iterate over times points
        
        # draw O[i,q,t] from mixture distribution (Eq. 3)
        indicator = sample(c(FALSE,TRUE),1,prob=c(1-gamma_val[[q]][t,which(model$prior.theta[i,]==1)],gamma_val[[q]][t,which(model$prior.theta[i,]==1)]))
        O[i,q,t]= ifelse(indicator,rdens_1(),abs(rdens_0()))
        #O[i,q,t]= indicator
      } 
    }
    
  }
  
  return(O)
}
#==========================================================

# Certaint degree of edge existence --------------------------- 
alpha_val = function(rho){
  
  return(matrix(rbern(length(rho),rho),nrow=nrow(rho),ncol = ncol(rho)
                #,dimnames =  list(netPara[["S_names"]],netPara[["S_names"]])
  ))       
}
##==========================================================

# edge-wise prior of W ------------------------------------
true.w= function(alpha,prior.W){
  return(  alpha*matrix(rlnorm(nrow(alpha)^2,meanlog = prior.W$nu,sdlog = prior.W$tau),nrow=nrow(alpha),ncol = ncol(alpha))+
             (1-alpha)*matrix(rtnorm(nrow(alpha)^2,mean = prior.W$nu0,sd = prior.W$tau0,lower = 0,upper = 0.2),nrow=nrow(alpha),ncol = ncol(alpha)))    #spike and slab prior 
}
##==========================================================
# motif rfunc (uVal as perturbation function)-------------------
rfunc <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    
    import <- uVal(t)
    #diag(parms[[1]])= diag(parms[[1]])*parms[[3]]
    dx= t(parms[[1]]) %*% (x) + (rep(import,ncol(parms[[1]])) * (t(parms[[4]])%*% parms[[2]]))
    
    list(dx)
  })
}
#######################################################################################################################
# Main
#######################################################################################################################

#initial parameters of network-----------------------------
S_num=4
P_num=4
E_num= 4
spSt=c("Bi","Dia","FF","FB","PC","JK","UVB")
tstep <- 0.1
#===========================================================


#constructing Model Object---------------------------------------------------------

PS<- diag(1,S_num) #P to S connections


l <- rep(list(0:1), S_num)
Q = as.matrix(expand.grid(l))    # known state of intervention nodes (interventing one Pnode at each experiment)
Q <- Q[-c(1,8),]

X=diag(x = 1,nrow =S_num )
prior.theta = do.call(rbind,lapply(1:S_num ,function(s)rbind(matrix(X[s,],nrow=E_num,ncol=S_num ,byrow=TRUE))))  #S_E connections

r0=rep(0,S_num)    #start state of hidden nodes
#r0 = c(0.2, 10, 2, 0, 0, 20, 0, 0, 0 , 20, 0) # UVB

prior.W=list(rho=NULL, nu0=NULL, tau0=NULL, nu=NULL, tau=NULL)  #weight matrix
#true_model$prior.W$lambda=1   #parameter for exponential function(nonedge)
prior.W$nu0=-0.3        #parameter for lnorm function(edge)lmean
prior.W$tau0=0.1    #parameter for lnorm function(edge)lsd
prior.W$nu=-0.2        #parameter for lnorm function(edge)lmean
prior.W$tau=0.5    #parameter for lnorm function(edge)lsd



true_model = list(W=NULL, Q=Q, prior.theta=prior.theta, prior.W=prior.W, orig.dat=NULL, r0=r0, PS=PS)



#sample a network Structure
# hidden_net=sampleSubnetwork(S_num)
graphSt= samplenetwork(S_num,P_num,E_num,true_model$prior.theta,stype = "Dia")
plot.igraph(graphSt) #plot graph of S, P and E nodes
#save(graphSt, file = paste0('graphSt.RData'))

true_model$prior.W$rho= as.matrix(get.adjacency(graphSt))[1:S_num,1:S_num] #it should be changed in case of prior knowledge of edges


alpha=alpha_val(true_model$prior.W$rho)

true_model$W=true.w(alpha,true_model$prior.W)
true_model$W[which(abs(true_model$W)<0.2,arr.ind = TRUE)]=0  #even small values for non edges can effect likelihood
diag(true_model$W)=0    # other situation 0, -0.2 because of gamma values 
#save(true_model, file = paste0('wt_true_model.RData'))



#===========================================================
#for(sim in 1:10)
#{
#simulate data ---------------------
true_model$orig.dat=sim.data(model = true_model, nE = E_num*S_num,times = 5, stype = "Dia", ts = tstep)

  
  #===========================================================
  
  #plot Observations-------------------
  
  my_palette <- colorRampPalette(c("blue","gray","red"))
  for (i in 1:nrow(Q)) {
    heatmap.2((true_model$orig.dat[,i,]),density.info = "none",col=my_palette,trace="none",Rowv = "none",Colv = "none",scale="none",main=paste("Exp_",i))
    
  }
  #d <- dev.off()
  #===========================================================
  
  true_ll<-Likelihood(true_model,true_model$orig.dat,compute_states ,gammas,ddens_1,ddens_0,true_model$PS,tstep)
  #save(l, file = paste0(sim,'_real_ll.RData'))
  
  #run MCMC and plot results------------------------------------
  res_mcmc <-learn(dat = true_model$orig.dat, Q = Q, rho=NULL, hyper.W=true_model$prior.W, prior.theta =true_model$prior.theta,seed=123,P_S = true_model$PS,ts = tstep)


plot(unlist(res_mcmc[[2]])) #all_ll
points(unlist(res_mcmc[[5]]),col="red") #accepted_ll
points(unlist(res_mcmc[[7]]),col="blue") #rejected_ll
points(true_ll,col="green") #true_ll


plot(x=c(1:length(res_mcmc[[5]]), y=unlist(res_mcmc[[5]])),type = "l") #trace of accepted_ll



plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][1,1]}),type = "l") #trace of acc w for each edge
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][1,2]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][1,3]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][2,1]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][2,2]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][2,3]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][3,1]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][3,2]}),type = "l")
plot(1:length(res_mcmc[[1]]),sapply(1:length(res_mcmc[[1]]), function(i){res_mcmc[[1]][[i]][3,3]}),type = "l")



postMatrices = postProcess(res_mcmc[[4]],true_model)


ROC_curve(postPrior = postMatrices[[1]],sc = postMatrices[[2]] ,model = true_model)


PR_curve(postMatrices[[1]],postMatrices[[2]] ,true_model)  
