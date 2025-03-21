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
  
  for (i in 1:nrow(mu)){

    sign=sign(mu[,1])


    mu[,1]<-mu[,1]*sign
    mu[,2]<-mu[,2]*sign
    mu[,3]<-mu[,3]*sign
  }
  
  
  return(mu)
}


# VMF_Regression_Old<-function(Data_list,iters=5000){
#   library(neurobase)
#   library(dti)
#   library(abind)
#   library(RNiftyReg)
#   library(rgl)
#   library(matlib)
#   library(igraph)
#   library(Rfast)
#   library(pracma)
#   library(gtools)
#   library(rotasym)
#   library(MASS)
#   library(Matrix)
#   library(ar.matrix)
#   
#   
#   #### Exact List Object to environment
#   list2env(Data_list, globalenv())
#   I=length(Y); V=length(Fiber_Ind); p=ncol(X); K=length(unique(Fiber_Ind))
#   
#   
#   #### Initial Values
#   sigma_Z<-sigma_Q<-1
#   kappa=rep(100,2)
#   rho<-rep(0.5,K)
#   omega<-rep(0.5,K)
#   alpha=matrix(NA,ncol=V,nrow=p); beta=matrix(NA,ncol=V,nrow=p)
#   gamma=matrix(0,ncol=V,nrow=p); nu=matrix(0,ncol=V,nrow=p)
#   
#   for (v in 1:V){
#     response=mu2theta_phi_scaled(t(sapply(1:I, function(i) Y[[i]][v,])))
#     alpha[,v]<-lm(response[,1]~0+X)$coef
#     beta[,v]<-lm(response[,2]~0+X)$coef
#   }
#   
#   
#   
#   
#   big_X<-lapply(1:I, function(i) bdiag(lapply(1:V, function(v) matrix(X[i,],nrow=1))))
#   theta_phi_scaled_all<-lapply(1:I, function(i) 
#     as.matrix(cbind(
#       big_X[[i]]%*%as.vector(alpha)+big_X[[i]]%*%as.vector(gamma)*(group_index[i]==1),
#       big_X[[i]]%*%as.vector(beta)+big_X[[i]]%*%as.vector(nu)*(group_index[i]==1)
#     ))
#   )
#   mu_all<-lapply(1:I, function(i) theta_phi_scaled2mu(theta_phi_scaled_all[[i]]))
#   Q_z=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k),sigma_Z,rho[k]) ))
#   Q_q=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k),sigma_Q,omega[k]) ))
#   
#   
#   
#   
#   kappa_MCMC=matrix(NA,iters,2)
#   alpha_MCMC=array(NA,dim=c(iters,p,V))
#   beta_MCMC=array(NA,dim=c(iters,p,V))
#   gamma_MCMC=array(NA,dim=c(iters,p,V))
#   nu_MCMC=array(NA,dim=c(iters,p,V))
#   rho_MCMC=matrix(NA,iters,K)
#   omega_MCMC=matrix(NA,iters,K)
#   sigma_Z_MCMC=matrix(NA,iters,1)
#   sigma_Q_MCMC=matrix(NA,iters,1)
#   
#   pb <- txtProgressBar(min = 0, max = iters, style = 3)
#   for (it in 1:iters){
#     setTxtProgressBar(pb, it)
#     
#     ##### Update kappa #####
#     for(group in 1:2){
#       
#       
#       kappa_old=exp(0+log(kappa[group]))
#       kappa_can=exp(rnorm(1,0,0.05)+log(kappa[group]))
#       like_old=sum(sapply(which(group_index==(group-1) ), function(i) sapply(1:V,function(v) d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(mu_all[[i]][v,]), kappa_old, log = TRUE))))
#       like_can=sum(sapply(which(group_index==(group-1) ), function(i) sapply(1:V,function(v) d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(mu_all[[i]][v,]), kappa_can, log = TRUE))))
#       
#       prob=min(exp(like_can-like_old),1)
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       
#       if(Accept==1){
#         kappa[group]=kappa_can
#       }
#     }
#     
#     kappa_MCMC[it,]<-kappa
#     
#     
#     ##### Update alpha #####
#     for(v in 1:V){
#       for(pp in 1:p){
#         
#         alpha_can<-alpha_old<-alpha
#         alpha_can[pp,v]<-alpha_can[pp,v]+rnorm(1,0,0.05)
#         
#         
#         old_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,])
#         can_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,]-c(X[i,]%*%alpha_old[,v]-X[i,]%*%alpha_can[,v],0) )
#         
#         old_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(old_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*alpha_old[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_z[[Fiber_Ind[v]]]%*%(alpha_old[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         can_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*alpha_can[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_z[[Fiber_Ind[v]]]%*%(alpha_can[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         prob=min(1,as.numeric(exp(can_like-old_like)))
#         Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#         if(Accept==1){
#           alpha<-alpha_can
#           for(i in 1:I){
#             theta_phi_scaled_all[[i]][v,]<-can_theta_phi_scaled_all[[i]]
#             mu_all[[i]][v,]<-theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])
#           }
#         }
#         
#         
#       }
#       
#     }
#     alpha_MCMC[it,,]<-alpha
#     
#     
#     
#     
#     
#     
#     ##### Update beta #####
#     for(v in 1:V){
#       for(pp in 1:p){
#         
#         beta_can<-beta_old<-beta
#         beta_can[pp,v]<-beta_can[pp,v]+rnorm(1,0,0.05)
#         
#         
#         old_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,])
#         can_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,]-c(0,X[i,]%*%beta_old[,v]-X[i,]%*%beta_can[,v]) )
#         
#         old_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(old_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*beta_old[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_q[[Fiber_Ind[v]]]%*%(beta_old[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         can_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*beta_can[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_q[[Fiber_Ind[v]]]%*%(beta_can[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         prob=min(1,as.numeric(exp(can_like-old_like)))
#         Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#         if(Accept==1){
#           beta<-beta_can
#           for(i in 1:I){
#             theta_phi_scaled_all[[i]][v,]<-can_theta_phi_scaled_all[[i]]
#             mu_all[[i]][v,]<-theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])
#           }
#         }
#         
#         
#       }
#       
#     }
#     beta_MCMC[it,,]<-beta
#     
#     
#     ##### Update gamma #####
#     for(v in 1:V){
#       for(pp in 1:p){
#         
#         gamma_can<-gamma_old<-gamma
#         gamma_can[pp,v]<-gamma_can[pp,v]+rnorm(1,0,0.05)
#         
#         
#         old_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,])
#         can_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,]-c(X[i,]%*%gamma_old[,v]-X[i,]%*%gamma_can[,v],0)*(group_index[i]==1) )
#         
#         old_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(old_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*gamma_old[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_z[[Fiber_Ind[v]]]%*%(gamma_old[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         can_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*gamma_can[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_z[[Fiber_Ind[v]]]%*%(gamma_can[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         prob=min(1,as.numeric(exp(can_like-old_like)))
#         Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#         if(Accept==1){
#           gamma<-gamma_can
#           for(i in 1:I){
#             theta_phi_scaled_all[[i]][v,]<-can_theta_phi_scaled_all[[i]]
#             mu_all[[i]][v,]<-theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])
#           }
#         }
#         
#         
#       }
#       
#     }
#     gamma_MCMC[it,,]<-gamma
#     
#     
#     
#     ##### Update nu #####
#     for(v in 1:V){
#       for(pp in 1:p){
#         
#         nu_can<-nu_old<-nu
#         nu_can[pp,v]<-nu_can[pp,v]+rnorm(1,0,0.05)
#         
#         
#         old_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,])
#         can_theta_phi_scaled_all=lapply(1:I,function(i) theta_phi_scaled_all[[i]][v,]-c(0,X[i,]%*%nu_old[,v]-X[i,]%*%nu_can[,v])*(group_index[i]==1) )
#         
#         old_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(old_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*nu_old[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_q[[Fiber_Ind[v]]]%*%(nu_old[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         can_like=sum(sapply(1:I, function(i)
#           d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])), kappa[group_index[i]+1], log = TRUE)))+
#           -0.5*nu_can[pp,Fiber_Ind==Fiber_Ind[v]]%*%Q_q[[Fiber_Ind[v]]]%*%(nu_can[pp,Fiber_Ind==Fiber_Ind[v]])
#         
#         
#         prob=min(1,as.numeric(exp(can_like-old_like)))
#         Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#         if(Accept==1){
#           nu<-nu_can
#           for(i in 1:I){
#             theta_phi_scaled_all[[i]][v,]<-can_theta_phi_scaled_all[[i]]
#             mu_all[[i]][v,]<-theta_phi_scaled2mu(can_theta_phi_scaled_all[[i]])
#           }
#         }
#         
#         
#       }
#       
#     }
#     nu_MCMC[it,,]<-nu
#     
#     
#     ##### Update rho #####
#     for (k in 1:K){
#       rho_old=rho[k]
#       rho_can=inv.logit(logit(rho_old)+rnorm(1,0,0.05))
#       Q_z_old=Q_z[[k]]
#       Q_z_can=Q.AR1(sum(Fiber_Ind==k),sigma_Z,rho_can)
#       
#       
#       like_old=sum(sapply(1:p, function(pp)   as.numeric(-0.5*alpha[pp,Fiber_Ind==k]%*%Q_z_old%*%alpha[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(-0.5*gamma[pp,Fiber_Ind==k]%*%Q_z_old%*%gamma[pp,Fiber_Ind==k])))+
#         -2*p*0.5*log(det(Q_z_old))
#       like_can=sum(sapply(1:p, function(pp)   as.numeric(-0.5*alpha[pp,Fiber_Ind==k]%*%Q_z_can%*%alpha[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(-0.5*gamma[pp,Fiber_Ind==k]%*%Q_z_can%*%gamma[pp,Fiber_Ind==k])))+
#         -2*p*0.5*log(det(Q_z_can))
#       
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         rho[k]=rho_can
#         Q_z[[k]]=Q_z_can
#       }
#     }
#     rho_MCMC[it,]<-rho
#     
#     
#     
#     
#     ##### Update omega #####
#     for (k in 1:K){
#       omega_old=omega[k]
#       omega_can=inv.logit(logit(omega_old)+rnorm(1,0,0.05))
#       Q_q_old=Q_q[[k]]
#       Q_q_can=Q.AR1(sum(Fiber_Ind==k),sigma_Q,omega_can)
#       
#       
#       like_old=sum(sapply(1:p, function(pp)   as.numeric(-0.5*beta[pp,Fiber_Ind==k]%*%Q_q_old%*%beta[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(-0.5*nu[pp,Fiber_Ind==k]%*%Q_q_old%*%nu[pp,Fiber_Ind==k])))+
#         -2*p*0.5*log(det(Q_q_old))
#       like_can=sum(sapply(1:p, function(pp)   as.numeric(-0.5*beta[pp,Fiber_Ind==k]%*%Q_q_can%*%beta[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(-0.5*nu[pp,Fiber_Ind==k]%*%Q_q_can%*%nu[pp,Fiber_Ind==k])))+
#         -2*p*0.5*log(det(Q_q_can))
#       
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         omega[k]=omega_can
#         Q_q[[k]]=Q_q_can
#       }
#     }
#     omega_MCMC[it,]<-omega
#     
#     
#     ##### sigma2_Z #######
#     a=sum(sapply(1:K, function(k)
#       sum(sapply(1:p, function(pp)   as.numeric(alpha[pp,Fiber_Ind==k]%*%(Q_z[[k]]*sigma_Z^2)%*%alpha[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(gamma[pp,Fiber_Ind==k]%*%(Q_z[[k]]*sigma_Z^2)%*%gamma[pp,Fiber_Ind==k])))
#     ))
#     b=2*p*V
#     sigma_Z<-sqrt(1/rgamma(1,a/2+5,b/2+5))
#     sigma_Z_MCMC[it,]<-sigma_Z
#     
#     
#     
#     
#     ##### sigma2_Q #######
#     a=sum(sapply(1:K, function(k)
#       sum(sapply(1:p, function(pp)   as.numeric(beta[pp,Fiber_Ind==k]%*%(Q_q[[k]]*sigma_Q^2)%*%beta[pp,Fiber_Ind==k])))+
#         sum(sapply(1:p, function(pp)   as.numeric(nu[pp,Fiber_Ind==k]%*%(Q_q[[k]]*sigma_Q^2)%*%nu[pp,Fiber_Ind==k])))
#     ))
#     b=2*p*V
#     sigma_Q<-sqrt(1/rgamma(1,a/2+5,b/2+5))
#     sigma_Q_MCMC[it,]<-sigma_Q
#     
#     
#     
#     
#     
#   }
#   
#   
#   return(list(kappa_MCMC=kappa_MCMC,
#               alpha_MCMC=alpha_MCMC,
#               beta_MCMC=beta_MCMC,
#               gamma_MCMC=gamma_MCMC,
#               nu_MCMC=nu_MCMC,
#               rho_MCMC=rho_MCMC,
#               omega_MCMC=omega_MCMC,
#               sigma_Z_MCMC=sigma_Z_MCMC,
#               sigma_Q_MCMC=sigma_Q_MCMC))
#   
# }
# 
# 
# 
# VMF_Regression_Old2<-function(Data_list,iters=5000){
#   library(neurobase)
#   library(dti)
#   library(abind)
#   library(RNiftyReg)
#   library(rgl)
#   library(matlib)
#   library(igraph)
#   library(Rfast)
#   library(pracma)
#   library(gtools)
#   library(rotasym)
#   library(MASS)
#   library(Matrix)
#   library(ar.matrix)
#   
#   
#   #### Exact List Object to environment
#   list2env(Data_list, globalenv())
#   I=length(Y); V=length(Fiber_Ind); p=ncol(X); K=length(unique(Fiber_Ind)); G=length(unique(group_index))
#   
#   
#   #### Initial Values
#   sigma_theta<-sigma_phi<-1
#   tau_theta<-tau_phi<-1
#   kappa=100
#   rho_alpha<-rep(0.1,K)
#   rho_beta<-rep(0.1,K)
#   rho_epsilon<-rep(0.1,K)
#   rho_xi<-rep(0.1,K)
#   
#   
#   
#   alpha=lapply(1:G, function(g) matrix(NA,ncol=V,nrow=p))
#   beta=lapply(1:G, function(g) matrix(NA,ncol=V,nrow=p))
#   
#   
#   
#   for (v in 1:V){
#     for(g in 1:G){
#       index=which(group_index==g)
#       response=mu2theta_phi_scaled(t(sapply(index, function(i) Y[[i]][v,])))
#       alpha[[g]][,v]<-lm(response[,1]~0+X[index,])$coef
#       beta[[g]][,v]<-lm(response[,2]~0+X[index,])$coef
#     }
#   }
#   
#   mu_alpha=lapply(1:G, function(g)  rowMeans(alpha[[g]]))
#   mu_beta=lapply(1:G, function(g)  rowMeans(beta[[g]]))
#   
#   
#   
#   big_X<-lapply(1:I, function(i) bdiag(lapply(1:V, function(v) matrix(X[i,],nrow=1))))
#   Mean_theta_phi_all<-lapply(1:I, function(i)  cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
#                                                big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
#   mu_all<-Y
#   theta_phi_scaled_all<-lapply(1:I, function(i) mu2theta_phi_scaled(mu_all[[i]]))
#   theta_phi_scaled_all_check<-theta_phi_scaled_all
#   
#   Q_resi_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_theta,rho_epsilon[k]) ))
#   Q_resi_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_phi,rho_xi[k]) ))
#   
#   Q_coef_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_theta,rho_alpha[k]) ))
#   Q_coef_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_phi,rho_beta[k]) ))
#   
#  
#   
#   #### MCMC Updates
#   kappa_MCMC=matrix(NA,iters,1)
#   alpha_MCMC=lapply(1:G,function(g) array(NA,dim=c(iters,p,V)) )
#   beta_MCMC=lapply(1:G,function(g) array(NA,dim=c(iters,p,V)) )
#  
#   rho_epsilon_MCMC=matrix(NA,iters,K)
#   rho_xi_MCMC=matrix(NA,iters,K)
#   rho_alpha_MCMC=matrix(NA,iters,K)
#   rho_beta_MCMC=matrix(NA,iters,K)
#   
#   
#   sigma_theta_MCMC=matrix(NA,iters,1)
#   sigma_phi_MCMC=matrix(NA,iters,1)
#   tau_theta_MCMC=matrix(NA,iters,1)
#   tau_phi_MCMC=matrix(NA,iters,1)
#   
#   theta_phi_scaled_all_MCMC=lapply(1:iters, function(it) lapply(1:I, function(i) NA*mu2theta_phi_scaled(mu_all[[i]]) ))
#   
#   
#   pb <- txtProgressBar(min = 0, max = iters, style = 3)
#   for (it in 1:iters){
#     setTxtProgressBar(pb, it)
#     
#     
#     ##### Metropolis Hastings:#######
#     ##### Update kappa:       #######
#     kappa_old=exp(0+log(kappa))
#     kappa_can=exp(rnorm(1)+log(kappa))
#     like_old=sum(sapply(1:I, function(i) sapply(1:V,function(v) d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(mu_all[[i]][v,]), kappa_old, log = TRUE))))
#     like_can=sum(sapply(1:I, function(i) sapply(1:V,function(v) d_vMF(matrix(Y[[i]][v,],nrow=1), as.vector(mu_all[[i]][v,]), kappa_can, log = TRUE))))
#     prob=min(exp(like_can-like_old),1)
#     Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#     if(Accept==1){
#       kappa<-kappa_can
#     }
# 
#     ### Store MCMC update
#     kappa_MCMC[it,]<-kappa
#     
#     
#     ##### Metropolis Hastings:#######
#     ##### theta_phi_scaled_all:#######
#     for (i in 1:I){
#       for (v in 1:V){
# 
#         theta_phi_scaled_old<-theta_phi_scaled_can<-theta_phi_scaled_all[[i]]
#         theta_phi_scaled_can[v,]<-theta_phi_scaled_old[v,]+rnorm(2)
# 
# 
#         Resi_theta_phi_old<-theta_phi_scaled_old-Mean_theta_phi_all[[i]]
#         Resi_theta_phi_can<-theta_phi_scaled_can-Mean_theta_phi_all[[i]]
# 
# 
# 
#         old_like=d_vMF(Y[[i]][v,], as.numeric(theta_phi_scaled2mu(theta_phi_scaled_old[v,])) , kappa,TRUE)+
#           -0.5*Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],1]%*%Q_resi_theta[[Fiber_Ind[v]]]%*%(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],1])+
#           -0.5*Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],2]%*%Q_resi_phi[[Fiber_Ind[v]]]%*%(Resi_theta_phi_old[Fiber_Ind==Fiber_Ind[v],2])
#         can_like=d_vMF(Y[[i]][v,], as.numeric(theta_phi_scaled2mu(theta_phi_scaled_can[v,])) , kappa,TRUE)+
#           -0.5*Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],1]%*%Q_resi_theta[[Fiber_Ind[v]]]%*%(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],1])+
#           -0.5*Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],2]%*%Q_resi_phi[[Fiber_Ind[v]]]%*%(Resi_theta_phi_can[Fiber_Ind==Fiber_Ind[v],2])
# 
#         prob=min(1,as.numeric(exp(can_like-old_like)))
#         Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#         if(Accept==1){
#             theta_phi_scaled_all[[i]]<-theta_phi_scaled_can
#             mu_all[[i]]<-theta_phi_scaled2mu(theta_phi_scaled_can)
#         }
# 
#       }
#     }
#     
#     
#     ### Store MCMC update
#     theta_phi_scaled_all_MCMC[[it]]<-theta_phi_scaled_all
#     
#     
#     ##### Gibbs Sampling:#######
#     ##### alpha:#######
#     for(g in 1:G){
#       index=which(group_index==g)
#       B=lapply(1:length(index), function(i) NULL)
#       b=lapply(1:length(index), function(i) NULL)
#       ii=0
#       for (i in index){
#          ii=ii+1
#          B[[ii]]=(t(big_X[[i]])%*%bdiag(Q_resi_theta)%*%big_X[[i]]+ bdiag(Q_coef_theta)%x%diag(p))
#          b[[ii]]=t(big_X[[i]])%*%bdiag(Q_resi_theta)%*%theta_phi_scaled_all[[i]][,1]+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_alpha[[g]] )
#       }
#       
#       
#       
# 
#       #+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_alpha[[g]] )
# 
#       MEAN=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
#       VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
#       # alpha[[g]]=matrix(t(chol(VAR))%*%(MEAN+rnorm(p*V)),nrow=p)
#       alpha[[g]]=matrix(rmvnorm(1,MEAN,VAR),nrow=p)
#       Mean_theta_phi_all<-lapply(1:I, function(i)  cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
#                                                          big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
# 
#       alpha_MCMC[[g]][it,,]<-alpha[[g]]
# 
# 
#     }
# 
# 
# 
# 
# 
#     ##### beta:#######
#     for(g in 1:G){
#       index=which(group_index==g)
#       B=lapply(1:length(index), function(i) NULL)
#       b=lapply(1:length(index), function(i) NULL)
#       ii=0
#       for (i in index){
#         ii=ii+1
#         B[[ii]]=(t(big_X[[i]])%*%bdiag(Q_resi_phi)%*%big_X[[i]]+ bdiag(Q_coef_phi)%x%diag(p))
#         b[[ii]]=t(big_X[[i]])%*%bdiag(Q_resi_phi)%*%theta_phi_scaled_all[[i]][,2]+(bdiag(Q_coef_phi)%x%diag(p))%*%matrix( rep(1,V) %x%mu_beta[[g]] )
#       }
# 
#       #+(bdiag(Q_coef_theta)%x%diag(p))%*%matrix( rep(1,V) %x%mu_beta[[g]] )
# 
#       MEAN=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))%*%Reduce("+",lapply(1:length(index),function(ii) (b[[ii]])))
#       VAR=solve(Reduce("+",lapply(1:length(index),function(ii) (B[[ii]]))))
#       # beta[[g]]=matrix(t(chol(VAR))%*%(MEAN+rnorm(p*V)),nrow=p)
#       beta[[g]]=matrix(rmvnorm(1,MEAN,VAR),nrow=p)
#       Mean_theta_phi_all<-lapply(1:I, function(i)  cbind(big_X[[i]]%*%as.vector(alpha[[group_index[i]]]),
#                                                          big_X[[i]]%*%as.vector(beta[[group_index[i]]])))
# 
#       beta_MCMC[[g]][it,,]<-beta[[g]]
#     }
#     
#     
#     
#     ##### Metropolis-Hastings:#######
#     ##### Update rho_alpha #####
#     for (k in 1:K){
#       rho_alpha_old=rho_alpha[k]
#       rho_alpha_can=inv.logit(logit(rho_alpha_old)+rnorm(1)*0.05)
#       Q_coef_theta_old=Q_coef_theta[[k]]
#       Q_coef_theta_can=Q.AR1(sum(Fiber_Ind==k),sigma_theta,rho_alpha_can)
# 
#       alpha_resi=lapply(1:G, function(g) alpha[[g]]-mu_alpha[[g]])
# 
# 
#       like_old=sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(-0.5*alpha_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_theta_old%*%alpha_resi[[g]][pp,Fiber_Ind==k])))))+
#         -G*p*0.5*log(det(Q_coef_theta_old))
#       like_can=sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp)   as.numeric(-0.5*alpha_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_theta_can%*%alpha_resi[[g]][pp,Fiber_Ind==k])))))+
#         -G*p*0.5*log(det(Q_coef_theta_can))
# 
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         rho_alpha[k]=rho_alpha_can
#         Q_coef_theta[[k]]=Q_coef_theta_can
#       }
#     }
#     rho_alpha_MCMC[it,]<-rho_alpha
# 
# 
#     ##### Metropolis-Hastings:#######
#     ##### Update rho_beta #####
#     for (k in 1:K){
#       rho_beta_old=rho_beta[k]
#       rho_beta_can=inv.logit(logit(rho_beta_old)+rnorm(1)*0.05)
#       Q_coef_phi_old=Q_coef_phi[[k]]
#       Q_coef_phi_can=Q.AR1(sum(Fiber_Ind==k),sigma_phi,rho_beta_can)
# 
#       beta_resi=lapply(1:G, function(g) beta[[g]]-mu_beta[[g]])
# 
# 
#       like_old=sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(-0.5*beta_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_phi_old%*%beta_resi[[g]][pp,Fiber_Ind==k])))))+
#         -G*p*0.5*log(det(Q_coef_phi_old))
#       like_can=sum(sapply(1:G,function(g) sum(sapply(1:p, function(pp)   as.numeric(-0.5*beta_resi[[g]][pp,Fiber_Ind==k]%*%Q_coef_phi_can%*%beta_resi[[g]][pp,Fiber_Ind==k])))))+
#         -G*p*0.5*log(det(Q_coef_phi_can))
# 
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         rho_beta[k]=rho_beta_can
#         Q_coef_phi[[k]]=Q_coef_phi_can
#       }
#     }
#     rho_beta_MCMC[it,]<-rho_beta
# 
#     ##### Gibbs Sampling:#######
#     ##### Update sigma_theta  sigma_phi #####
#     alpha_resi=lapply(1:G, function(g) alpha[[g]]-mu_alpha[[g]])
#     a=sum(sapply(1:K, function(k) sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(alpha_resi[[g]][pp,Fiber_Ind==k]%*%(Q_coef_theta[[k]]*sigma_theta^2)%*%alpha_resi[[g]][pp,Fiber_Ind==k]))))))
#     b=G*p*V
#     sigma_theta<-sqrt(1/rgamma(1,a/2+1,b/2+1))
#     
#     beta_resi=lapply(1:G, function(g) beta[[g]]-mu_beta[[g]])
#     a=sum(sapply(1:K, function(k) sapply(1:G,function(g) sum(sapply(1:p, function(pp) as.numeric(beta_resi[[g]][pp,Fiber_Ind==k]%*%(Q_coef_phi[[k]]*sigma_phi^2)%*%beta_resi[[g]][pp,Fiber_Ind==k]))))))
#     b=G*p*V
#     sigma_phi<-sqrt(1/rgamma(1,a/2+1,b/2+1))
#     
#     Q_coef_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_theta,rho_alpha[k]) ))
#     Q_coef_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), sigma_phi,rho_beta[k]) ))
#     
# 
#     ##### Metropolis-Hastings:#######
#     ##### Update rho_epsilon  rho_xi #####
#     Resi_theta_phi_all<-lapply(1:I,function(i) theta_phi_scaled_all[[i]]-Mean_theta_phi_all[[i]])
#     Resi_theta=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,1])
#     Resi_phi=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,2])
# 
#     for (k in 1:K){
#       rho_epsilon_old=rho_epsilon[k]
#       rho_epsilon_can=inv.logit(logit(rho_epsilon_old)+rnorm(1)*0.05)
#       Q_resi_theta_old=Q_resi_theta[[k]]
#       Q_resi_theta_can=Q.AR1(sum(Fiber_Ind==k),tau_theta,rho_epsilon_can)
# 
# 
# 
#       like_old=sum(sapply(1:I,function(i) as.numeric(-0.5*Resi_theta[Fiber_Ind==k,i]%*%Q_resi_theta_old%*%Resi_theta[Fiber_Ind==k,i])))+
#         -I*0.5*log(det(Q_resi_theta_old))
#       like_can=sum(sapply(1:I,function(i) as.numeric(-0.5*Resi_theta[Fiber_Ind==k,i]%*%Q_resi_theta_can%*%Resi_theta[Fiber_Ind==k,i])))+
#         -I*0.5*log(det(Q_resi_theta_can))
# 
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         rho_epsilon[k]=rho_epsilon_can
#         Q_resi_theta[[k]]=Q_resi_theta_can
#       }
#     }
#     rho_epsilon_MCMC[it,]<-rho_epsilon
# 
# 
#     for (k in 1:K){
#       rho_xi_old=rho_xi[k]
#       rho_xi_can=inv.logit(logit(rho_xi_old)+rnorm(1)*0.05)
#       Q_resi_phi_old=Q_resi_phi[[k]]
#       Q_resi_phi_can=Q.AR1(sum(Fiber_Ind==k),tau_phi,rho_xi_can)
# 
# 
# 
#       like_old=sum(sapply(1:I,function(i) as.numeric(-0.5*Resi_phi[Fiber_Ind==k,i]%*%Q_resi_phi_old%*%Resi_phi[Fiber_Ind==k,i])))+
#         -I*0.5*log(det(Q_resi_phi_old))
#       like_can=sum(sapply(1:I,function(i) as.numeric(-0.5*Resi_phi[Fiber_Ind==k,i]%*%Q_resi_phi_can%*%Resi_phi[Fiber_Ind==k,i])))+
#         -I*0.5*log(det(Q_resi_phi_can))
# 
#       prob=min(1,exp(like_can-like_old))
#       Accept=sample(c(0,1),size=1,prob=c(1-prob,prob))
#       if(Accept==1){
#         rho_xi[k]=rho_xi_can
#         Q_resi_phi[[k]]=Q_resi_phi_can
#       }
#     }
#     rho_xi_MCMC[it,]<-rho_xi
# 
# 
# 
#     ##### Gibbs Sampling:#######
#     ##### Update tau_theta  tau_phi #####
#     Resi_theta_phi_all<-lapply(1:I,function(i) theta_phi_scaled_all[[i]]-Mean_theta_phi_all[[i]])
#     Resi_theta=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,1])
#     Resi_phi=sapply(1:I, function(i) Resi_theta_phi_all[[i]][,2])
# 
#     a=sum(sapply(1:K,function(k) sapply(1:I,function(i) as.numeric(Resi_theta[Fiber_Ind==k,i]%*% (Q_resi_theta[[k]]*tau_theta^2) %*%Resi_theta[Fiber_Ind==k,i]))))
#     b=I*V
#     tau_theta<-sqrt(1/rgamma(1,a/2+1,b/2+1))
# 
#     a=sum(sapply(1:K,function(k) sapply(1:I,function(i) as.numeric(Resi_phi[Fiber_Ind==k,i]%*% (Q_resi_phi[[k]]*tau_phi^2) %*%Resi_phi[Fiber_Ind==k,i]))))
#     b=I*V
#     tau_phi<-sqrt(1/rgamma(1,a/2+1,b/2+1))
# 
#     Q_resi_theta=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_theta,rho_epsilon[k]) ))
#     Q_resi_phi=(lapply(1:K, function(k) Q.AR1(sum(Fiber_Ind==k), tau_phi,rho_xi[k]) ))
#     tau_theta_MCMC[it,]<-tau_theta
#     tau_phi_MCMC[it,]<-tau_phi
#     
#     
#   } 
#     
#     
#     
#     
#     
#     
#  
#   
#   return(list(kappa_MCMC=kappa_MCMC,
#               alpha_MCMC=alpha_MCMC,
#               beta_MCMC=beta_MCMC,
#               gamma_MCMC=gamma_MCMC,
#               nu_MCMC=nu_MCMC,
#               rho_MCMC=rho_MCMC,
#               omega_MCMC=omega_MCMC,
#               sigma_Z_MCMC=sigma_Z_MCMC,
#               sigma_Q_MCMC=sigma_Q_MCMC))
#   
# }
# 
# 
# 































