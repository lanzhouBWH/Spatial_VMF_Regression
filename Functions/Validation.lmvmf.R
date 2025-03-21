Validation.lmvmf<-function(fit,Data_list,iters=5000,burnins=0){
 
  dvmf<-function (y, mu, k, logden = FALSE) 
  {
    y <- as.matrix(y)
    p <- dim(y)[2]
    if (p == 1) 
      y <- t(y)
    p <- dim(y)[2]
    den <- (p/2 - 1) * log(k) - 0.5 * p * log(2 * pi) + k * tcrossprod(mu, 
                                                                       y) - log(besselI(k, p/2 - 1, expon.scaled = TRUE)) - 
      k
    if (!logden) 
      den <- exp(den)
    den
  }
  
  
  #### Exact List Object to environment#################################################################
  list2env(Data_list, globalenv())
  I=length(Y); V=length(Fiber_Ind); p=ncol(X); K=length(unique(Fiber_Ind)); G=length(unique(group_index))
  ######################################################################################################
  
  
  
  pred_Matrix=array(NA, dim=c(3,I,V,iters-burnins))
  
  index=0
  for(it in (burnins+1):iters ){
    
    index=index+1
    
    ### Exact values on iteration i
    ALPHA=lapply(1:G, function(g) fit$alpha_MCMC[[g]][it,,])
    BETA=lapply(1:G, function(g) fit$beta_MCMC[[g]][it,,])
    kappa=fit$kappa_MCMC[it]
    
    
    ### Simulate Residuals
    rho_epsilon=fit$rho_epsilon_MCMC[it,]
    rho_xi=fit$rho_xi_MCMC[it,]
    tau_theta=fit$tau_theta_MCMC[it,]
    tau_phi=fit$tau_phi_MCMC[it,]
    
    
    S_resi_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon,corr=TRUE)*tau_theta^2)    ))
    S_resi_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi,corr=TRUE)*tau_phi^2) ))
    
    EPSILON=sapply(1:I, function(i)
      unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_theta[[k]]))))
    XI=sapply(1:I, function(i)
      unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_phi[[k]]))))
    
    
    
    for(i in 1:I){
      for(v in 1:V){
        
        theta_scaled=X[i,]%*%ALPHA[[group_index[i]]][,v]+EPSILON[v,i]
        phi_scaled=X[i,]%*%BETA[[group_index[i]]][,v]+XI[v,i]
        
        #theta_scaled=X[i,]%*%ALPHA[[group_index[i]]][,v]
        #phi_scaled=X[i,]%*%BETA[[group_index[i]]][,v]
        mu=theta_phi_scaled2mu(c(theta_scaled,phi_scaled,1))
        
        pred_Matrix[,i,v,index]<-solve(fit$U_MCMC[[it]])%*%t(mu)
        #print(DIC_Matrix[i,v,index])
      }
    }
    
    
  }
  
  PP_Matrix=array(NA, dim=c(3,I,V))
  
  for(i in 1:I){
    PP_Matrix[,i,]<-sapply(1:V, function(v) vmf.mle(t(pred_Matrix[,i,v,]))$mu ) 
  }
  
  
  
  alpha=sapply(1:I, function(i) sapply(1:V, function(v) acos(PP_Matrix[,i,v]%*%Y[[i]][v,])))
  kappa=mean(fit$kappa_MCMC[(burnins+1):iters])
  Log_Lik=sapply(1:I, function(i) sapply(1:V, function(v) dvmf(Y[[i]][v,],as.numeric(PP_Matrix[,i,v]),kappa,logden = TRUE)))   
  
  
  RR_VMF=list(alpha=alpha,Log_Lik=Log_Lik )
  return(RR_VMF)
  
}