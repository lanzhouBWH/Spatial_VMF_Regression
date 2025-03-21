MCMC_Sampler<-function(data,X){
  
  niter=5000
  adaptive=5000
  kappa_MCMC=rep(NA,niter)
  alpha_MCMC=matrix(NA,2,niter)
  beta_MCMC=matrix(NA,2,niter)
  theta_phi_scaled_MCMC=array(NA,dim=c(nrow(data),2,niter))
  U_MCMC<-lapply(1:niter, function(it) matrix(NA,3,3)  )
  
  #### Defiining Initial Values #####
  kappa=Directional::vmf.mle(data)$kappa
  
  EST<-function(data,X){
    
    response=mu2theta_phi_scaled(data)
    alpha=coef(lm(response[,1]~0+X))
    alpha[is.na(alpha)]=0
    beta=coef(lm(response[,2]~0+X))
    beta[is.na(beta)]=0
    
    eee=theta_phi_scaled2mu(cbind(X%*%alpha,X%*%beta))
    
    kappa=Rfast::vmf.mle(data)$kappa
    return(list(alpha=alpha,beta=beta,kappa=kappa))
  }
  
  RR=EST(data,X)
  alpha=RR$alpha
  beta=RR$beta
  sigma_alpha=1; sigma_beta=1
  theta_phi_scaled=mu2theta_phi_scaled(data)

  
  skewvec <- rnorm(3)*0
  A <- matrix(0, 3, 3)
  A[upper.tri(A)] <- skewvec
  A <- A - t(A)
  U   <- (diag(3) + A) %*% solve(diag(3) - A)
  
  
  
  
  
  S_kappa=1
  S_alpha=rep(1,ncol(X))
  S_beta=rep(1,ncol(X))
  
  I=nrow(data)
  
  S_theta_phi=array(0,dim=c(I,2,2))
  for(i in 1:I){
      S_theta_phi[i,,]=diag(2)
  }
  
  S_skewvec=diag(3)
  
  for(it in 1:niter){
    
    
    #### U change #####
    u=S_skewvec%*%rnorm(3)
    skewvec_can=skewvec+u
    A <- matrix(0, 3, 3)
    A[upper.tri(A)] <- skewvec_can
    A <- A - t(A)
    U_can   <- t( (diag(3) - A) %*% solve(diag(3) + A) )
    
    theta_phi_scaled_can=mu2theta_phi_scaled(t(U_can%*%t(theta_phi_scaled2mu(theta_phi_scaled))))
    
    vvv_can<-sum(sapply(1:I,function(i) dnorm(theta_phi_scaled_can[i,1],(X%*%alpha)[i],sqrt(sigma_alpha),TRUE)))+
    sum(sapply(1:I,function(i) dnorm(theta_phi_scaled_can[i,2],(X%*%alpha)[i],sqrt(sigma_beta),TRUE)))
    
    vvv_old<-sum(sapply(1:I,function(i) dnorm(theta_phi_scaled[i,1],(X%*%alpha)[i],sqrt(sigma_alpha),TRUE)))+
      sum(sapply(1:I,function(i) dnorm(theta_phi_scaled[i,2],(X%*%alpha)[i],sqrt(sigma_beta),TRUE)))
    
    prob=min(1,exp(vvv_can-vvv_old))
    if(prob>runif(1)){
      theta_phi_scaled<-theta_phi_scaled_can
      U=U_can
    }
    
    if(it<=adaptive){
      S_skewvec <- ramcmc::adapt_S(S_skewvec, u, prob, it)
    }
    
    
    
    ### alpha ###
    post_alpha <- Gibbs.regression(X[,2], theta_phi_scaled[,1], NULL, 1, trace='bsmt', fix='xy',intercept=TRUE,
                                   Ngauss=0)
    alpha=as.numeric(post_alpha$B)
    sigma_alpha=as.numeric(post_alpha$Sigma)
    
    ### beta ###
    post_beta <- Gibbs.regression(X[,2], theta_phi_scaled[,2], NULL, 1, trace='bsmt', fix='xy',intercept=TRUE,
                                   Ngauss=0)
    beta=as.numeric(post_beta$B)
    sigma_beta=as.numeric(post_beta$Sigma)
    
    
    ### Theta
    for (i in 1:I){

      #my_k=Fiber_Ind[v]
      #my_v=V_Ktype[v]

      u=S_theta_phi[i,,]%*%rnorm(2)
      theta_phi_scaled_old<-theta_phi_scaled[i,]
      theta_phi_scaled_can<-theta_phi_scaled[i,]+ as.numeric(u)

      vvv_can=dvmf(data[i,], k=kappa, theta_phi_scaled2mu(theta_phi_scaled_can),logden = TRUE)+
        dnorm(theta_phi_scaled_can[1],(X%*%alpha)[i], sqrt(sigma_alpha),TRUE)+
        dnorm(theta_phi_scaled_can[2],(X%*%beta)[i], sqrt(sigma_beta),TRUE)
      vvv_old=dvmf(data[i,], k=kappa, theta_phi_scaled2mu(theta_phi_scaled_old),logden = TRUE)+
        dnorm(theta_phi_scaled_old[1],(X%*%alpha)[i], sqrt(sigma_alpha),TRUE)+
        dnorm(theta_phi_scaled_old[2],(X%*%beta)[i], sqrt(sigma_beta),TRUE)

      prob=min(1,exp(vvv_can-vvv_old))
      if(prob>runif(1)){
        theta_phi_scaled[i,]<-theta_phi_scaled_can
      }

      if(it<=adaptive){
        S_theta_phi[i,,] <- ramcmc::adapt_S(S_theta_phi[i,,], u, prob, it)
      }

    }
  
    
    ##### Metropolis Hasting for Kappa ###########
    kappa_old=exp(0+log(kappa))
    u=S_kappa*rnorm(1)
    kappa_can=max(0.1,exp(u+log(kappa)))
    eee=theta_phi_scaled2mu(theta_phi_scaled)
    
    vvv_can=sum(sapply(1:nrow(data),function(i) dvmf(data[i,], k=kappa_can, eee[i,],logden = TRUE)))
    vvv=sum(sapply(1:nrow(data),function(i) dvmf(data[i,], k=kappa, eee[i,],logden = TRUE)))
    if(exp(vvv_can-vvv)>runif(1)){
      kappa=kappa_can
    }
    
    if (it<=adaptive){
      S_kappa<-ramcmc::adapt_S(S_kappa, u, min(1,exp(vvv_can-vvv)), it)
    }
    
    kappa_MCMC[it]<-kappa
    alpha_MCMC[,it]<-alpha
    beta_MCMC[,it]<-beta
    theta_phi_scaled_MCMC[,,it]<-theta_phi_scaled
    U_MCMC[[it]]<-U
    

  }
  
  
  
  return(list(U_MCMC=U_MCMC,
              kappa_MCMC=kappa_MCMC,
              alpha_MCMC=alpha_MCMC,
              beta_MCMC=beta_MCMC,
              theta_phi_scaled_MCMC=theta_phi_scaled_MCMC))
  
  
}



Off2Direction<-function(C){
  Matrix=diag(3)
  Matrix[1,1]<-Matrix[1,1]<-C[1]
  Matrix[2,1]<-Matrix[1,2]<-C[2]
  Matrix[3,1]<-Matrix[1,3]<-C[3]
  Matrix[2,2]<-Matrix[2,2]<-C[4]
  Matrix[3,2]<-Matrix[3,2]<-C[5]
  Matrix[3,3]<-Matrix[3,3]<-C[5]
  
  DATA=as.numeric(eigen(Matrix)$vectors[,1])
  sign=sign(DATA[1])
  
  
  DATA[1]<-DATA[1]*sign
  DATA[2]<-DATA[2]*sign
  DATA[3]<-DATA[3]*sign
  
  return(DATA)
}



Cond_Prod1<-function(Resi,v,Fiber_Ind,rho,s,V){
  
  if(v==1){
    result=dnorm(Resi[v],rho*Resi[v+1],sqrt(1-rho^2)*s,log = TRUE)
  }else if(Fiber_Ind[v-1]!=Fiber_Ind[v]){
    result=dnorm(Resi[v],rho*Resi[v+1],sqrt(1-rho^2)*s,log = TRUE)
  }else if(v==V){
    result=dnorm(Resi[v-1],rho*Resi[v],sqrt(1-rho^2)*s,log = TRUE)+
      dnorm(Resi[v],rho*Resi[v],s,log = TRUE)
  }else if(Fiber_Ind[v]!=Fiber_Ind[v+1]){
    result=dnorm(Resi[v-1],rho*Resi[v],sqrt(1-rho^2)*s,log = TRUE)+
      dnorm(Resi[v],rho*Resi[v],s,log = TRUE)
  }else{
    result=dnorm(Resi[v],rho*Resi[v+1],sqrt(1-rho^2)*s,log = TRUE)+
      dnorm(Resi[v-1],rho*Resi[v],sqrt(1-rho^2)*s,log = TRUE)
  }
  return(result)
  
}


mu2theta_phi_scaled<-function(mu){
  
  if(!is.matrix(mu)){mu=matrix(mu,nrow=1)}
  theta_phi<-cart2sph(mu)
  if(!is.matrix(theta_phi)){theta_phi=matrix(theta_phi,nrow=1)}
  theta_phi<-theta_phi[,1:2]
  if(!is.matrix(theta_phi)){theta_phi=matrix(theta_phi,nrow=1)}
  return(cbind(logit((theta_phi[,1]+pi)/(2*pi)),logit((theta_phi[,2]+pi/2)/(pi))))
}



theta_phi_scaled2mu<-function(theta_phi_scaled){
  
  if(!is.matrix(theta_phi_scaled)){theta_phi_scaled=matrix(theta_phi_scaled,nrow=1)}
  theta_phi<-cbind(inv.logit(theta_phi_scaled[,1])*(2*pi)-pi,inv.logit(theta_phi_scaled[,2])*(pi)-pi/2)
  mu<-sph2cart(cbind(theta_phi,1))
  
  
  
  if(!is.matrix(mu)){mu=matrix(mu,nrow=1)}

  # for (i in 1:nrow(mu)){
  # 
  #   sign=sign(mu[,1])
  # 
  # 
  #   mu[,1]<-mu[,1]*sign
  #   mu[,2]<-mu[,2]*sign
  #   mu[,3]<-mu[,3]*sign
  # }
  
  
  return(mu)
}
