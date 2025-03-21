Synthetic_Data<-function(kappa,X,group_index,COEF,Fiber_Ind,P){
  
  I=nrow(X); V=ncol(COEF$alpha[[1]]); p=ncol(X);G=length(unique(group_index));K=length(unique(Fiber_Ind));
  
  rho_alpha<-inla.ar.pacf2phi(runif(P))
  rho_beta<-inla.ar.pacf2phi(runif(P))
  rho_epsilon<-inla.ar.pacf2phi(runif(P))
  rho_xi<-inla.ar.pacf2phi(runif(P))
  
  S_resi_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon,corr=TRUE)*0.5^2)    ))
  S_resi_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi,corr=TRUE)*0.5^2) ))
  
  S_coef_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_alpha,corr=TRUE)*0.5^2) ))
  S_coef_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_beta,corr=TRUE)*0.5^2) ))
  
  ALPHA=lapply(1:G, function(g) COEF$alpha[[g]])
  BETA=lapply(1:G, function(g) COEF$beta[[g]])
  
  
  Y=list()
  
  for(i in 1:I){
    
    EPSILON=unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_theta[[k]])))
    XI=unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_phi[[k]])))
    
    temp=NULL
    
    for(v in 1:V){
      
      theta_scaled=X[i,]%*%ALPHA[[group_index[i]]][,v]+EPSILON[v]
      phi_scaled=X[i,]%*%BETA[[group_index[i]]][,v]+XI[v]
      
      mu=rvmf(1,as.numeric(theta_phi_scaled2mu(c(theta_scaled,phi_scaled,1))),kappa)
      temp=rbind(temp, mu)
      
      
    }
    
    Y[[i]]<-temp
    
    
  }
  
  Data_list=list(Y=Y, Fiber_Ind=Fiber_Ind,group_index=group_index, X=X)
  
  return(Data_list)
  
  
}