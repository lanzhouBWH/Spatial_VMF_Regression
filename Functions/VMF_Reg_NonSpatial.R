VMF_Reg_NonSpatial<-function(Data_list,iters=50000,burnins=0,adaptive=200,P=3,fit=NULL,save.dir=NULL,Saving_Freq=50,Para=TRUE){
 
  
  
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
  
  
  #### Register Cluster ################################################################################
  #myCluster <- makeCluster(detectCores()-1, # number of cores to use
  #                         type = "SOCK") # type of cluster
  #registerDoParallel(myCluster)
  if(Para==TRUE){
    registerDoParallel(cores=detectCores()-1)
  }
  ######################################################################################################
  
  
  if(!is.null(fit)){
    
    list2env(fit, globalenv())
    
  }else{
    
    it=1
    
    #### Exact List Object to environment#################################################################
    list2env(Data_list, globalenv())
    I=length(Y); V=length(Fiber_Ind); p=ncol(X); K=length(unique(Fiber_Ind)); G=length(unique(group_index))
    ######################################################################################################
    
    
    #### Conditional Distribution Sepecifcation ##########################################################
    #COND_DIST=lapply(1:K, function(k)  vecchia_specify(matrix(1:sum(Fiber_Ind==k),ncol=1),m=P)$U.prep$revNNarray)
    #V_Ktype=unlist(lapply(1:K, function(k) 1:sum(Fiber_Ind==k) ))
    
    #### Initial Values (Will be updated) ################################################################
    
    ## Kappa ######################################
    kappa<-vmf.mle(list.rbind(Y))$kappa
    ###############################################
    
    ## Correlation Parameters #####################
    rho_alpha<-inla.ar.pacf2phi(runif(P)*0)
    rho_beta<-inla.ar.pacf2phi(runif(P)*0)
    rho_epsilon<-inla.ar.pacf2phi(runif(P)*0)
    rho_xi<-inla.ar.pacf2phi(runif(P)*0)
    ###############################################
    
    
    ## Coefficients Parameters ############################
    alpha=lapply(1:G, function(g) matrix(NA,ncol=V,nrow=p))
    beta=lapply(1:G, function(g) matrix(NA,ncol=V,nrow=p))
    vv_theta_temp=NULL
    vv_phi_temp=NULL
    for (v in 1:V){
      for(g in 1:G){
        index=which(group_index==g)
        response=mu2theta_phi_scaled(t(sapply(index, function(i) Y[[i]][v,])))
        lmobj1=lm(response[,1]~0+X[index,])
        lmobj2=lm(response[,2]~0+X[index,])
        alpha[[g]][,v]<-lmobj1$coef
        beta[[g]][,v]<-lmobj2$coef
        
        
        ### if NA, then 0 #############
        alpha[[g]][,v][is.na(alpha[[g]][,v])]=0
        beta[[g]][,v][is.na(beta[[g]][,v])]=0
        ##############################
        
        
        vv_theta_temp=c(vv_theta_temp,summary(lmobj1)$sigma)
        vv_phi_temp=c(vv_theta_temp,summary(lmobj2)$sigma)
      }
    }
    mu_alpha=lapply(1:G, function(g)  rowMeans(alpha[[g]]))
    mu_beta=lapply(1:G, function(g)  rowMeans(beta[[g]]))
    #######################################################
    
    
    ## Initial Value of Variances #########################
    sigma_theta<-tau_theta<-1
    sigma_phi<-tau_phi<-1
    #######################################################
    
    
    ## Initial Values of Transformed Means ################
    big_X<-lapply(1:I, function(i) diag(V)%x%matrix(X[i,],nrow=1) )
    Mean_theta_phi_all<-lapply(1:I, function(i)  as.matrix(cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
                                                                 big_X[[i]]%*%as.vector(beta[[group_index[i]]])))  )
    mu_all<-Y
    theta_phi_scaled_all<-lapply(1:I, function(i) mu2theta_phi_scaled(mu_all[[i]]))
    theta_phi_scaled_all_check<-theta_phi_scaled_all
    
    mu_all<-lapply(1:I, function(i) theta_phi_scaled2mu(cbind(Mean_theta_phi_all[[i]],1)))
    #######################################################
    
    
    
    ## Initial Values of Transformed Means ################
    S_resi_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon,corr=TRUE)*tau_theta^2)    ))
    S_resi_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi,corr=TRUE)*tau_phi^2) ))
    
    Q_resi_theta=(lapply(1:K, function(k) solve(S_resi_theta[[k]])    ))
    Q_resi_phi=(lapply(1:K, function(k) solve(S_resi_phi[[k]]) ))
    
    S_coef_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_alpha,corr=TRUE)*sigma_theta^2) ))
    S_coef_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_beta,corr=TRUE)*sigma_phi^2) ))
    
    Q_coef_theta=(lapply(1:K, function(k) solve(S_coef_theta[[k]])    ))
    Q_coef_phi=(lapply(1:K, function(k) solve(S_coef_phi[[k]]) ))
    
    #######################################################
    
    
    ##### Initial Values for 
    skewvec <- rnorm(3)*0
    A <- matrix(0, 3, 3)
    A[upper.tri(A)] <- skewvec
    A <- A - t(A)
    U   <- (diag(3) + A) %*% solve(diag(3) - A)
    
    
    
    #############################################################################################################
    
    
    #### MCMC Updates Storing ####################################################################
    kappa_MCMC=matrix(NA,iters-burnins,1)
    alpha_MCMC=lapply(1:G,function(g) array(NA,dim=c(iters-burnins,p,V)) )
    beta_MCMC=lapply(1:G,function(g) array(NA,dim=c(iters-burnins,p,V)) )
    
    rho_epsilon_MCMC=matrix(NA,iters-burnins,P)
    rho_xi_MCMC=matrix(NA,iters-burnins,P)
    rho_alpha_MCMC=matrix(NA,iters-burnins,P)
    rho_beta_MCMC=matrix(NA,iters-burnins,P)
    
    sigma_theta_MCMC=matrix(NA,iters-burnins,1)
    sigma_phi_MCMC=matrix(NA,iters-burnins,1)
    tau_theta_MCMC=matrix(NA,iters-burnins,1)
    tau_phi_MCMC=matrix(NA,iters-burnins,1)
    
    #theta_phi_scaled_all_MCMC=lapply(1:(iters-burnins), function(it) lapply(1:I, function(i) NA*mu2theta_phi_scaled(mu_all[[i]]) ))
    
    U_MCMC<-lapply(1:(iters-burnins), function(it) matrix(NA,3,3)  )
    ################################################################################################
    
    
    ######## Step Size for Metropolis Hastings #####################################################
    S_kappa=0.1
    M_kappa=0.1
    
    S_rho_alpha=rep(1,P)*0.1
    S_rho_beta=rep(1,P)*0.1
    S_rho_epsilon=rep(1,P)*0.1
    S_rho_xi=rep(1,P)*0.1
    
    
    S_skewvec=diag(3)*0.1
    
    
    S_theta_phi=array(0,dim=c(I,V,2,2))
    for(i in 1:I){
      for(v in 1:V){
        S_theta_phi[i,v,,]=1*diag(2)
      }
    }
    M=array(1,dim=c(I,V,2))
    
    ###############################################################################################
  }
  
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  for (it in it:iters){
    setTxtProgressBar(pb, it)
    
    ##### Metropolis Hastings:#######
    ##### theta_phi_scaled_all:#######
    OOO<-foreach(i=1:I ) %dopar% {
      
      ii=sample(1:I,1)
      
      if(TRUE){
        
        
        for (v in 1:V){
          
          #my_k=Fiber_Ind[v]
          #my_v=V_Ktype[v]
          
          u=S_theta_phi[i,v,,]%*%rnorm(2)
          theta_phi_scaled_old<-theta_phi_scaled_can<-theta_phi_scaled_all[[i]]
          theta_phi_scaled_can[v,]<-theta_phi_scaled_old[v,]+ u
          
          
          Resi_theta_phi_old<-theta_phi_scaled_old-Mean_theta_phi_all[[i]]
          Resi_theta_phi_can<-theta_phi_scaled_can-Mean_theta_phi_all[[i]]
          
          like_old=dvmf(U%*%Y[[i]][v,],  as.numeric(theta_phi_scaled2mu(theta_phi_scaled_old[v,]),kappa), TRUE)+
            as.numeric(-0.5*(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],1]%*%Q_resi_theta[[Fiber_Ind[v]]])%*%(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],1]))+
            as.numeric(-0.5*(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],2]%*%Q_resi_phi[[Fiber_Ind[v]]])%*%(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],2]))
          like_can=dvmf(U%*%Y[[i]][v,],as.numeric(theta_phi_scaled2mu(theta_phi_scaled_can[v,]),kappa),TRUE)+
            as.numeric(-0.5*(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],1]%*%Q_resi_theta[[Fiber_Ind[v]]])%*%(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],1]))+
            as.numeric(-0.5*(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],2]%*%Q_resi_phi[[Fiber_Ind[v]]])%*%(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],2]))
          
          # like_old=dvmf(Y[[i]][v,],  as.numeric(theta_phi_scaled2mu(theta_phi_scaled_old[v,]),kappa),TRUE)
          # like_can=dvmf(Y[[i]][v,],as.numeric(theta_phi_scaled2mu(theta_phi_scaled_can[v,]),kappa),TRUE)
          # 
          # x1=Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],1]
          # y1=Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],1]
          # x2=Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],2]
          # y2=Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],2]
          # A=Q_resi_theta[[Fiber_Ind[v]]]
          
          
          prob=min(1,as.numeric(exp(like_can-like_old
          )))
          # Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
          MH=(like_can-like_old)
          Accept<-(log(runif(1))<MH)
          if (is.na(Accept)|is.na(prob)){Accept=FALSE; prob=0}
          if(Accept){
            
            theta_phi_scaled_all[[i]]<-theta_phi_scaled_can
            mu_all[[i]][v,]<-theta_phi_scaled2mu(theta_phi_scaled_can[v,])
          }
          
          if(it<=adaptive){
            S_theta_phi[i,v,,] <- ramcmc::adapt_S(S_theta_phi[i,v,,], u, prob, it)
          }
          
        }
      }
      out=list(theta_phi_scaled_all[[i]],mu_all[[i]],S_theta_phi[i,,,])
    }
    
    theta_phi_scaled_all=lapply(1:I, function(i) OOO[[i]][[1]])
    mu_all=lapply(1:I, function(i) OOO[[i]][[2]])
    for (i in 1:I){
      S_theta_phi[i,,,]<-OOO[[i]][[3]]
    }
    
    
    if(it>burnins){
      ### Store MCMC update
      #theta_phi_scaled_all_MCMC[[it-burnins]]<-theta_phi_scaled_all
    } 
    #############################################################################################
    
    
    
    #### U change #####
    u=S_skewvec%*%rnorm(3)
    skewvec_can=skewvec+u
    A <- matrix(0, 3, 3)
    A[upper.tri(A)] <- skewvec_can
    A <- A - t(A)
    U_can   <-  t( (diag(3) - A) %*% solve(diag(3) + A))
    
    like_old=0
    like_can=0
    for(i in 1:I){
      
      
      for(v in 1:V){
        
        ttt=theta_phi_scaled2mu(theta_phi_scaled_all[[i]][v,])
        like_old=like_old+dvmf(U%*%Y[[i]][v,],  ttt,kappa, TRUE)
        like_can=like_can+dvmf(U_can%*%Y[[i]][v,],  ttt,kappa, TRUE)
        
        
        
      }
    }
    
    prob=min(1,exp(like_can-like_old))
    if(prob>runif(1)){
      U=U_can
      skewvec=skewvec_can
      #print(U)
    }
    
    if(it<=adaptive){
      S_skewvec <- ramcmc::adapt_S(S_skewvec, u, prob, it)
    }
    
    if(it>burnins){
      ### Store MCMC update
      U_MCMC[[it-burnins]]<-U
    }
    ########################################################################################
    
    
    
    
    
    
    ##### Metropolis-Hastings:###################################################################
    ##### Update rho_alpha ######################################################################
    alpha_resi=lapply(1:G, function(g) alpha[[g]]-mu_alpha[[g]])
    for (PP in 1:P){
      rho_alpha_old<-rho_alpha
      pacf_alpha_old<-inla.ar.phi2pacf(rho_alpha_old)
      #u<-rnorm(1)*S_rho_alpha[PP]
      pacf_alpha_can<-pacf_alpha_old
      #pacf_alpha_can[PP]<-min(0.99,inv.logit(logit(pacf_alpha_can[PP])+u))
      pacf_alpha_can[PP]<-tanh(rnorm(1,0,0.5))*0
      rho_alpha_can<-inla.ar.pacf2phi(pacf_alpha_can)
      
      
      Q_coef_theta_old=Q_coef_theta
      Q_coef_theta_can<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_alpha_can,corr=TRUE)*sigma_theta^2) ) )
      
      
      #like_old<-dbeta(pacf_alpha_old[PP],10,1,log=TRUE)
      #like_can<-dbeta(pacf_alpha_can[PP],10,1,log=TRUE)
      
      like_old<-0
      like_can<-0
      for (k in 1:K){
        like_old=like_old+sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) 0.5*log(det(Q_coef_theta_old[[k]]))+as.numeric(-0.5*alpha_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_theta_old[[k]]%*%alpha_resi[[g]][pp,Fiber_Ind==k])))))
        like_can=like_can+sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) 0.5*log(det(Q_coef_theta_can[[k]]))+as.numeric(-0.5*alpha_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_theta_can[[k]]%*%alpha_resi[[g]][pp,Fiber_Ind==k])))))
      }
      
      
      prob=min(1,exp(like_can-like_old))
      # Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
      MH=(like_can-like_old)
      Accept<-(log(runif(1))<MH)
      if (is.na(Accept)|is.na(prob)){Accept=FALSE;prob=0}
      if(Accept){
        rho_alpha=rho_alpha_can
        Q_coef_theta=Q_coef_theta_can
      }
      
      # if(it<=adaptive){
      #   S_rho_alpha[PP] <- ramcmc::adapt_S(S_rho_alpha[PP], u, prob, it+1)
      # }
      if (it>burnins){
        rho_alpha_MCMC[it-burnins,]<-rho_alpha
      }
      
    }
    ####################################################################################
    
    ##### Metropolis-Hastings:##########################################################
    ##### Update rho_beta ##############################################################
    beta_resi=lapply(1:G, function(g) beta[[g]]-mu_beta[[g]])
    for (PP in 1:P){
      rho_beta_old<-rho_beta
      pacf_beta_old<-inla.ar.phi2pacf(rho_beta_old)
      #u<-rnorm(1)*S_rho_beta[PP]
      pacf_beta_can<-pacf_beta_old
      #pacf_beta_can[PP]<-min(0.99,inv.logit(logit(pacf_beta_can[PP])+u))
      pacf_beta_can[PP]<-tanh(rnorm(1,0,0.5))*0
      rho_beta_can<-inla.ar.pacf2phi(pacf_beta_can)
      
      Q_coef_phi_old=Q_coef_phi
      Q_coef_phi_can<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_beta_can,corr=TRUE)*sigma_phi^2) ) )
      
      
      #like_old<-dbeta(pacf_beta_old[PP],10,1,log=TRUE)
      #like_can<-dbeta(pacf_beta_can[PP],10,1,log=TRUE)
      
      like_old<-0
      like_can<-0
      
      for (k in 1:K){
        like_old=like_old+sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) 0.5*log(det(Q_coef_phi_old[[k]]))+as.numeric(-0.5*beta_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_phi_old[[k]]%*%beta_resi[[g]][pp,Fiber_Ind==k])))))
        like_can=like_can+sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) 0.5*log(det(Q_coef_phi_can[[k]]))+as.numeric(-0.5*beta_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_phi_can[[k]]%*%beta_resi[[g]][pp,Fiber_Ind==k])))))
      }
      
      
      prob=min(1,exp(like_can-like_old))
      # Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
      MH=(like_can-like_old)
      Accept<-(log(runif(1))<MH)
      if (is.na(Accept)|is.na(prob)){Accept=FALSE; prob=0}
      if(Accept){
        rho_beta=rho_beta_can
        Q_coef_phi=Q_coef_phi_can
      }
      
      # if(it<=adaptive){
      #   S_rho_beta[PP] <- ramcmc::adapt_S(S_rho_beta[PP], u, prob, it+1)
      # }
      if (it>burnins){
        rho_beta_MCMC[it-burnins,]<-rho_beta
      }
      
    }
    ########################################################################################
    
    
    ##### Gibbs Sampling:###################################################################
    ##### Update sigma_theta  sigma_phi ####################################################
    alpha_resi=lapply(1:G, function(g) alpha[[g]]-mu_alpha[[g]])
    a=sum(sapply(1:K, function(k) sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(alpha_resi[[g]][pp,Fiber_Ind==k]%*%(Q_coef_theta[[k]]*sigma_theta^2)%*%alpha_resi[[g]][pp,Fiber_Ind==k]))))))
    b=G*p*V
    sigma_theta<-sqrt(1/rgamma(1,a/2+0.1,b/2+0.1))
    
    beta_resi=lapply(1:G, function(g) beta[[g]]-mu_beta[[g]])
    a=sum(sapply(1:K, function(k) sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(beta_resi[[g]][pp,Fiber_Ind==k]%*%(Q_coef_phi[[k]]*sigma_phi^2)%*%beta_resi[[g]][pp,Fiber_Ind==k]))))))
    b=G*p*V
    sigma_phi<-sqrt(1/rgamma(1,a/2+0.1,b/2+0.1))
    
    # Q_coef_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_theta,rho_alpha) ))
    # Q_coef_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_phi,rho_beta) ))
    
    Q_coef_theta<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_alpha,corr=TRUE)*sigma_theta^2) ))
    Q_coef_phi<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_beta,corr=TRUE)*sigma_phi^2) ))
    
    
    if(it>burnins){
      sigma_theta_MCMC[it-burnins,]<-sigma_theta
      sigma_phi_MCMC[it-burnins,]<-sigma_phi
    }
    #############################################################################################
    
    
    ##### Metropolis-Hastings:###################################################################
    ##### Update rho_epsilon  rho_xi ############################################################
    ##### PART0: Preparing Residuals ############################################################
    Resi_theta_phi_all<-lapply(1:I,function(i) theta_phi_scaled_all[[i]]-Mean_theta_phi_all[[i]])
    Resi_theta=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,1])
    Resi_phi=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,2])
    #############################################################################################
    ##### PART1.1. Updating rho_epsilon #########################################################
    for (PP in 1:P){
      rho_epsilon_old<-rho_epsilon
      pacf_epsilon_old<-inla.ar.phi2pacf(rho_epsilon_old)
      #u<-rnorm(1)*S_rho_epsilon[PP]
      pacf_epsilon_can<-pacf_epsilon_old
      #pacf_epsilon_can[PP]<-min(0.99,inv.logit(logit(pacf_epsilon_can[PP])+u))
      pacf_epsilon_can[PP]<-tanh(rnorm(1,0,0.5))*0
      rho_epsilon_can<-inla.ar.pacf2phi(pacf_epsilon_can)
      
      
      Q_resi_theta_old=Q_resi_theta
      Q_resi_theta_can<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon_can,corr=TRUE)*tau_theta^2)    ))
      
      like_old<-0
      like_can<-0
      for(k in 1:K){
        like_old=like_old+sum(sapply(1:I,function(i) 0.5*log(det(Q_resi_theta_old[[k]]))+as.numeric(-0.5*Resi_theta[Fiber_Ind==k,i]%*%Q_resi_theta_old[[k]]%*%Resi_theta[Fiber_Ind==k,i])))
        like_can=like_can+sum(sapply(1:I,function(i) 0.5*log(det(Q_resi_theta_can[[k]]))+as.numeric(-0.5*Resi_theta[Fiber_Ind==k,i]%*%Q_resi_theta_can[[k]]%*%Resi_theta[Fiber_Ind==k,i])))
      }
      
      
      prob=min(1,exp(like_can-like_old));prob
      # Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
      MH=(like_can-like_old)
      Accept<-(log(runif(1))<MH)
      if (is.na(Accept)|is.na(prob)){Accept=FALSE;prob=0}
      if(Accept){
        rho_epsilon=rho_epsilon_can
        Q_resi_theta=Q_resi_theta_can
      }
      
      # if(it<=adaptive){
      #   S_rho_epsilon[PP] <- ramcmc::adapt_S(S_rho_epsilon[PP], u, prob, it+1)
      # }
      if(it>burnins){
        rho_epsilon_MCMC[it-burnins,]<-rho_epsilon
      }
      
      
    }
    ##############################################################################################
    ##### PART1.2. Updating rho_xi ###############################################################
    for (PP in 1:P){
      rho_xi_old<-rho_xi
      pacf_xi_old<-inla.ar.phi2pacf(rho_xi_old)
      #u<-rnorm(1)*S_rho_xi[PP]
      pacf_xi_can<-pacf_xi_old
      #pacf_xi_can[PP]<-min(0.99,inv.logit(logit(pacf_xi_can[PP])+u))
      pacf_xi_can[PP]<-tanh(rnorm(1,0,0.5))*0
      rho_xi_can<-inla.ar.pacf2phi(pacf_xi_can)
      
      Q_resi_phi_old=Q_resi_phi
      Q_resi_phi_can<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi_can,corr=TRUE)*tau_theta^2)    ))
      
      #like_old<-dbeta(pacf_xi_old[PP],10,1,log=TRUE  )
      #like_can<-dbeta(pacf_xi_can[PP],10,1,log=TRUE  )
      
      like_old<-0
      like_can<-0
      
      for(k in 1:K){
        like_old=like_old+sum(sapply(1:I,function(i) 0.5*log(det(Q_resi_phi_old[[k]]))+as.numeric(-0.5*Resi_phi[Fiber_Ind==k,i]%*%Q_resi_phi_old[[k]]%*%Resi_phi[Fiber_Ind==k,i])))
        like_can=like_can+sum(sapply(1:I,function(i) 0.5*log(det(Q_resi_phi_can[[k]]))+as.numeric(-0.5*Resi_phi[Fiber_Ind==k,i]%*%Q_resi_phi_can[[k]]%*%Resi_phi[Fiber_Ind==k,i])))
      }
      
      
      prob=min(1,exp(like_can-like_old));prob
      # Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
      MH=(like_can-like_old)
      Accept<-(log(runif(1))<MH)
      if (is.na(Accept)|is.na(prob)){Accept=FALSE;prob=0}
      if(Accept){
        rho_xi=rho_xi_can
        Q_resi_phi=Q_resi_phi_can
      }
      
      # if(it<=adaptive){
      #   S_rho_xi[PP] <- ramcmc::adapt_S(S_rho_xi[PP], u, prob, it+1)
      # }
      if(it>burnins){
        rho_xi_MCMC[it-burnins,]<-rho_xi
      }
      
    }
    #############################################################################################
    
    
    ##### Gibbs Sampling:#################################################################################
    ##### Update tau_theta  tau_phi ######################################################################
    ##### PART0: Preparing Residuals #####################################################################
    Resi_theta_phi_all<-lapply(1:I,function(i) theta_phi_scaled_all[[i]]-Mean_theta_phi_all[[i]])
    Resi_theta=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,1])
    Resi_phi=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,2])
    ######################################################################################################
    
    ##### PART1: Preparing Parameters #####################################################################
    a=sum(sapply(1:K,function(k) sapply(1:I,function(i) as.numeric(Resi_theta[Fiber_Ind==k,i]%*% (Q_resi_theta[[k]]*tau_theta^2) %*%Resi_theta[Fiber_Ind==k,i]))))
    b=I*V
    tau_theta<-sqrt(1/rgamma(1,a/2+0.1,b/2+0.1))
    
    a=sum(sapply(1:K,function(k) sapply(1:I,function(i) as.numeric(Resi_phi[Fiber_Ind==k,i]%*% (Q_resi_phi[[k]]*tau_phi^2) %*%Resi_phi[Fiber_Ind==k,i]))))
    b=I*V
    tau_phi<-sqrt(1/rgamma(1,a/2+0.1,b/2+0.1))
    
    #Q_resi_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_theta,rho_epsilon) ))
    #Q_resi_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_phi,rho_xi) ))
    
    Q_resi_theta<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon,corr=TRUE)*tau_theta^2)    ))
    Q_resi_phi<-(lapply(1:K, function(k) solve(ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi,corr=TRUE)*tau_phi^2) ))
    
    if(it>burnins){
      tau_theta_MCMC[it-burnins,]<-tau_theta
      tau_phi_MCMC[it-burnins,]<-tau_phi
    }
    ######################################################################################################
    
    
    
    
    ##### Metropolis Hastings:###################################################################
    ##### Update kappa:       ###################################################################
    kappa_old=exp(0+log(kappa))
    u=S_kappa*rnorm(1)
    kappa_can=max(20,exp(u+log(kappa)))
    #M_kappa_old=M_kappa
    #M_kappa_can=M_kappa_old+rnorm(1)
    like_old=sum(sapply(1:I, function(i) sapply(1:V,function(v) dvmf(U%*%Y[[i]][v,],  mu_all[[i]][v,], kappa_old, logden = TRUE))))
    like_can=sum(sapply(1:I, function(i) sapply(1:V,function(v) dvmf(U%*%Y[[i]][v,],  mu_all[[i]][v,], kappa_can, logden = TRUE))))
    prob=min(exp(like_can-like_old),1)
    #Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
    MH=(like_can-like_old)
    Accept<-(log(runif(1))<MH)
    if (is.na(Accept)|is.na(prob) ){Accept=FALSE; prob=0 }
    if(Accept){
      kappa<-kappa_can
      #M_kappa=-M_kappa_can
    }
    
    if (it<=adaptive){
      S_kappa<-ramcmc::adapt_S(S_kappa, u, prob, it)
    }
    if(it>burnins){
      ### Store MCMC update
      kappa_MCMC[it-burnins,]<-kappa
    }
    ###############################################################################################
    
    
    ##### Gibbs Sampling:##########################################################################
    ##### alpha:###################################################################################
    for(g in 1:G){
      for(k in 1:K){
        
        V_k=dim(Q_resi_theta[[k]])[1]
        V_index=which(Fiber_Ind==k)
        
        index=which(group_index==g)
        B=lapply(1:length(index), function(i) NULL)
        b=lapply(1:length(index), function(i) NULL)
        ii=0
        for (i in index){
          ii=ii+1
          X_k=big_X[[i]][1:V_k,1:(V_k*p)]
          B[[ii]]=(t(X_k)%*%Q_resi_theta[[k]]%*%X_k+ Q_coef_theta[[k]]%x%diag(p))
          b[[ii]]=t(X_k)%*%Q_resi_theta[[k]]%*%theta_phi_scaled_all[[i]][V_index,1]+(Q_coef_theta[[k]]%x%diag(p))%*%matrix( rep(1,V_k) %x%mu_alpha[[g]] )
        }
        
        VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
        MEAN=VAR%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
        
        alpha[[g]][,V_index]=matrix(Rfast::rmvnorm(1,MEAN,VAR),nrow=p)
      }
      
      
      
      
      if(it>burnins){
        alpha_MCMC[[g]][it-burnins,,]<-alpha[[g]]
      }
      
    }
    
    
    Mean_theta_phi_all<-lapply(1:I, function(i)  as.matrix(cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
                                                                 big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
    )
    
    
    # for(g in 1:G){
    #   index=which(group_index==g)
    #   B=lapply(1:length(index), function(i) NULL)
    #   b=lapply(1:length(index), function(i) NULL)
    #   ii=0
    #   for (i in index){
    #     ii=ii+1
    #     B[[ii]]=(t(big_X[[i]])%*%bdiag(Q_resi_theta)%*%big_X[[i]]+ bdiag(Q_coef_theta)%x%diag(p))
    #     b[[ii]]=t(big_X[[i]])%*%bdiag(Q_resi_theta)%*%theta_phi_scaled_all[[i]][,1]+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_alpha[[g]] )
    #   }
    #   
    #   
    #   
    #   
    #   #+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_alpha[[g]] )
    #   
    #   #MEAN=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
    #   VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
    #   MEAN=VAR%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
    #   # alpha[[g]]=matrix(t(chol(VAR))%*%(MEAN+rnorm(p*V)),nrow=p)
    #   alpha[[g]]=matrix(rmvnorm(1,MEAN,VAR),nrow=p)
    #   Mean_theta_phi_all<-lapply(1:I, function(i)  as.matrix(cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
    #                                                                big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
    #   )
    #   
    #   if(it>burnins){
    #     alpha_MCMC[[g]][it-burnins,,]<-alpha[[g]]
    #   }
    #   
    # }
    #########################################################################################################
    
    ##### Gibbs Sampling:###################################################################################
    ##### beta:#############################################################################################
    for(g in 1:G){
      for(k in 1:K){
        
        V_k=dim(Q_resi_phi[[k]])[1]
        V_index=which(Fiber_Ind==k)
        
        index=which(group_index==g)
        B=lapply(1:length(index), function(i) NULL)
        b=lapply(1:length(index), function(i) NULL)
        ii=0
        for (i in index){
          ii=ii+1
          X_k=big_X[[i]][1:V_k,1:(V_k*p)]
          B[[ii]]=(t(X_k)%*%Q_resi_phi[[k]]%*%X_k+ Q_coef_theta[[k]]%x%diag(p))
          b[[ii]]=t(X_k)%*%Q_resi_phi[[k]]%*%theta_phi_scaled_all[[i]][V_index,2]+(Q_coef_theta[[k]]%x%diag(p))%*%matrix( rep(1,V_k) %x%mu_alpha[[g]] )
        }
        
        VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
        MEAN=VAR%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
        
        beta[[g]][,V_index]=matrix(Rfast::rmvnorm(1,MEAN,VAR),nrow=p)
      }
      
      
      
      
      if(it>burnins){
        beta_MCMC[[g]][it-burnins,,]<-beta[[g]]
      }
      
    }
    
    
    Mean_theta_phi_all<-lapply(1:I, function(i)  as.matrix(cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
                                                                 big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
    )
    
    
    
    
    
    
    
    # for(g in 1:G){
    #   index=which(group_index==g)
    #   B=lapply(1:length(index), function(i) NULL)
    #   b=lapply(1:length(index), function(i) NULL)
    #   ii=0
    #   for (i in index){
    #     ii=ii+1
    #     B[[ii]]=(t(big_X[[i]])%*%bdiag(Q_resi_phi)%*%big_X[[i]]+ bdiag(Q_coef_phi)%x%diag(p))
    #     b[[ii]]=t(big_X[[i]])%*%bdiag(Q_resi_phi)%*%theta_phi_scaled_all[[i]][,2]+(bdiag(Q_coef_phi)%x%diag(p))%*%matrix( rep(1,V) %x%mu_beta[[g]] )
    #   }
    #   
    #   #+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_beta[[g]] )
    #   
    #   #MEAN=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
    #   VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
    #   MEAN=VAR%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
    #   # beta[[g]]=matrix(t(chol(VAR))%*%(MEAN+rnorm(p*V)),nrow=p)
    #   beta[[g]]=matrix(rmvnorm(1,MEAN,VAR),nrow=p)
    #   Mean_theta_phi_all<-lapply(1:I, function(i)  as.matrix(cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
    #                                                                big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
    #   )
    #   if(it>burnins){
    #     beta_MCMC[[g]][it-burnins,,]<-beta[[g]]
    #   }
    # }
    #################################################################################################
    
    
    
    ###### Update mu_alpha and mu_beta ##############################################################
    
    mu_alpha=lapply(1:G, function(g)  rowMeans(alpha[[g]]))
    mu_beta=lapply(1:G, function(g)  rowMeans(beta[[g]]))
    
    
    #################################################################################################
    
    
    
    
    
    
    
    
    if ( !is.null(save.dir) ){
      if (it%%Saving_Freq==0){
        fit=list(
          it=it,
          Y=Y, Fiber_Ind=Fiber_Ind,group_index=group_index, X=X,
          I=I, V=length(Fiber_Ind), p=ncol(X), K=length(unique(Fiber_Ind)), G=length(unique(group_index)),
          ## Kappa ######################################
          kappa=kappa,U=U,skewvec=skewvec,
          ###############################################
          
          ## Correlation Parameters #####################
          rho_alpha=rho_alpha,
          rho_beta=rho_beta,
          rho_epsilon=rho_epsilon,
          rho_xi=rho_xi,
          ###############################################
          
          ## Coefficients Parameters ############################
          alpha=alpha,
          beta=beta,
          mu_alpha=mu_alpha,
          mu_beta=mu_beta,
          #######################################################
          
          ## Initial Value of Variances #########################
          sigma_theta=sigma_theta,
          sigma_phi=sigma_phi,
          tau_theta=tau_theta,
          tau_phi=tau_phi,
          #######################################################
          
          
          ## Initial Values of Transformed Means ################
          big_X=big_X,
          Mean_theta_phi_all=Mean_theta_phi_all,
          theta_phi_scaled_all=theta_phi_scaled_all,
          mu_all=mu_all,
          #######################################################
          
          ## Initial Values of Transformed Means ################
          Q_resi_theta=Q_resi_theta,
          Q_resi_phi=Q_resi_phi,
          
          Q_coef_theta=Q_coef_theta,
          Q_coef_phi=Q_coef_phi,
          #######################################################
          
          
          U_MCMC=U_MCMC,
          kappa_MCMC=kappa_MCMC,
          alpha_MCMC=alpha_MCMC,
          beta_MCMC=beta_MCMC,
          rho_epsilon_MCMC=rho_epsilon_MCMC,
          rho_xi_MCMC=rho_xi_MCMC,
          rho_alpha_MCMC=rho_alpha_MCMC,
          rho_beta_MCMC=rho_beta_MCMC,
          sigma_theta_MCMC=sigma_theta_MCMC,
          sigma_phi_MCMC=sigma_phi_MCMC,
          tau_theta_MCMC=tau_theta_MCMC,
          tau_phi_MCMC=tau_phi_MCMC,
          
          
          
          
          ######## Step Size for Metropolis Hastings #####################################################
          S_skewvec=S_skewvec,
          S_kappa=S_kappa,
          M_kappa=M_kappa,
          S_rho_alpha=S_rho_alpha,
          S_rho_beta=S_rho_beta,
          S_rho_epsilon=S_rho_epsilon,
          S_rho_xi=S_rho_xi,
          S_theta_phi=S_theta_phi,
          M=M
          ###############################################################################################
          #,
          #theta_phi_scaled_all_MCMC=theta_phi_scaled_all_MCMC
        )
        save(fit,file = save.dir)
      }
    }
    
    
    
    
  } 
  
  
  
  fit=list(
    it=it,
    Y=Y, Fiber_Ind=Fiber_Ind,group_index=group_index, X=X,
    I=I, V=length(Fiber_Ind), p=ncol(X), K=length(unique(Fiber_Ind)), G=length(unique(group_index)),
    ## Kappa ######################################
    kappa=kappa,U=U,skewvec=skewvec,
    ###############################################
    
    ## Correlation Parameters #####################
    rho_alpha=rho_alpha,
    rho_beta=rho_beta,
    rho_epsilon=rho_epsilon,
    rho_xi=rho_xi,
    ###############################################
    
    ## Coefficients Parameters ############################
    alpha=alpha,
    beta=beta,
    mu_alpha=mu_alpha,
    mu_beta=mu_beta,
    #######################################################
    
    ## Initial Value of Variances #########################
    sigma_theta=sigma_theta,
    sigma_phi=sigma_phi,
    tau_theta=tau_theta,
    tau_phi=tau_phi,
    #######################################################
    
    
    ## Initial Values of Transformed Means ################
    big_X=big_X,
    Mean_theta_phi_all=Mean_theta_phi_all,
    theta_phi_scaled_all=theta_phi_scaled_all,
    mu_all=mu_all,
    #######################################################
    
    ## Initial Values of Transformed Means ################
    Q_resi_theta=Q_resi_theta,
    Q_resi_phi=Q_resi_phi,
    
    Q_coef_theta=Q_coef_theta,
    Q_coef_phi=Q_coef_phi,
    #######################################################
    
    
    U_MCMC=U_MCMC,
    kappa_MCMC=kappa_MCMC,
    alpha_MCMC=alpha_MCMC,
    beta_MCMC=beta_MCMC,
    rho_epsilon_MCMC=rho_epsilon_MCMC,
    rho_xi_MCMC=rho_xi_MCMC,
    rho_alpha_MCMC=rho_alpha_MCMC,
    rho_beta_MCMC=rho_beta_MCMC,
    sigma_theta_MCMC=sigma_theta_MCMC,
    sigma_phi_MCMC=sigma_phi_MCMC,
    tau_theta_MCMC=tau_theta_MCMC,
    tau_phi_MCMC=tau_phi_MCMC,
    
    
    
    
    ######## Step Size for Metropolis Hastings #####################################################
    S_skewvec=S_skewvec,
    S_kappa=S_kappa,
    M_kappa=M_kappa,
    S_rho_alpha=S_rho_alpha,
    S_rho_beta=S_rho_beta,
    S_rho_epsilon=S_rho_epsilon,
    S_rho_xi=S_rho_xi,
    S_theta_phi=S_theta_phi,
    M=M
    ###############################################################################################
    #,
    #theta_phi_scaled_all_MCMC=theta_phi_scaled_all_MCMC
  )
  
  #stopCluster(myCluster)
  stopImplicitCluster()
  
  return(fit)
}
