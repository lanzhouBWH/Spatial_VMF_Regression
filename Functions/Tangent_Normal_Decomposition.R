Tangent_Normal_Decomposition<-function(fit,coord,Fiber_Ind,X_ref,X_adj,gg,iters=5000,burnins=2000,brain_templete,deciles=9,
                                       image_folder,zoom=0.15,gif=TRUE
){
  
 
  V=length(Fiber_Ind); K=length(unique(Fiber_Ind));
  
  Ref_Posterior_Matrix=array(NA, dim=c(3,V,iters-burnins))
  Adj_Posterior_Matrix=array(NA, dim=c(3,V,iters-burnins))
  R_Posterior_Matrix=array(NA, dim=c(3,V,iters-burnins))
  t_Posterior_Matrix=array(NA, dim=c(V,iters-burnins))
  
  
  index=0
  for(it in (burnins+1):iters ){
    
    index=index+1
    
    ### Exact values on iteration i
    ALPHA=fit$alpha_MCMC[[gg]][it,,]
    BETA=fit$beta_MCMC[[gg]][it,,]
    kappa=fit$kappa_MCMC[it]
    
    
    ### Simulate Residuals
    rho_epsilon=fit$rho_epsilon_MCMC[it,]
    rho_xi=fit$rho_xi_MCMC[it,]
    tau_theta=fit$tau_theta_MCMC[it,]
    tau_phi=fit$tau_phi_MCMC[it,]
    
    
    S_resi_theta=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_epsilon,corr=TRUE)*tau_theta^2)    ))
    S_resi_phi=(lapply(1:K, function(k) (ARMA.var(n=sum(Fiber_Ind==k), ar=rho_xi,corr=TRUE)*tau_phi^2) ))
    
    EPSILON=unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_theta[[k]])))
    XI=unlist(lapply(1:K, function(k) mvrnorm(n = 1, mu=rep(0,sum(Fiber_Ind==k)), Sigma=S_resi_phi[[k]])))
    
    
    
    
    for(v in 1:V){
      
      theta_scaled=X_ref%*%ALPHA[,v]+EPSILON[v]
      phi_scaled=X_ref%*%BETA[,v]+XI[v]
      
      #theta_scaled=X_pred%*%ALPHA[[gg]][,v]
      #phi_scaled=X_pred%*%BETA[[gg]][,v]
      mu=theta_phi_scaled2mu(c(theta_scaled,phi_scaled,1))
      
      
      Ref_Posterior_Matrix[,v,index]<-mu
      
      
      
      theta_scaled=X_adj%*%ALPHA[,v]+EPSILON[v]
      phi_scaled=X_adj%*%BETA[,v]+XI[v]
      
      #theta_scaled=X_pred%*%ALPHA[[gg]][,v]
      #phi_scaled=X_pred%*%BETA[[gg]][,v]
      mu=theta_phi_scaled2mu(c(theta_scaled,phi_scaled,1))
      
      
      Adj_Posterior_Matrix[,v,index]<-mu
      
      
      Costheta=Adj_Posterior_Matrix[,v,index]%*%Ref_Posterior_Matrix[,v,index]
      t=as.numeric(sqrt(1-Costheta^2))
      
      
      R_Posterior_Matrix[,v,index]<-(Adj_Posterior_Matrix[,v,index]-Ref_Posterior_Matrix[,v,index]*sqrt(1-t^2))/t
      
      t_Posterior_Matrix[v,index]<-t
    }
    
    
  }
  
  
  T_est=rowMeans(t_Posterior_Matrix)
  Ref_est=sapply(1:V, function(v) rowMeans(Ref_Posterior_Matrix[,v,])/sqrt(sum((rowMeans(Ref_Posterior_Matrix[,v,]))^2))  )
  Adj_est=sapply(1:V, function(v) rowMeans(Adj_Posterior_Matrix[,v,])/sqrt(sum((rowMeans(Adj_Posterior_Matrix[,v,]))^2))  )
  R_est=sapply(1:V, function(v) rowMeans(R_Posterior_Matrix[,v,])/sqrt(sum((rowMeans(R_Posterior_Matrix[,v,]))^2))  )
  
  
  
  dd=quantile(T_est, probs = seq(.1, .9, by = .1))
  selected=which(T_est>dd[deciles])
  
  
  
  
  ###################################################################################
  ###Begin: To Plot m values ########################################################
  ###################################################################################
  # open3d()
  # # Needs to be a bigger window than the default
  # par3d(windowRect = c(100*5, 100*5, 612*5, 612*5))
  # 
  # 
  # 
  # 
  # 
  # for(v in 1:(V-1)){
  #   if(Fiber_Ind[v]==Fiber_Ind[v+1]){
  #     rgl.linestrips(c(coord[v,1], coord[v+1,1]),
  #                    c(coord[v,2], coord[v+1,2]),
  #                    c(coord[v,3], coord[v+1,3]),col="grey")
  #   }
  # }
  # 
  # cols <- myColorRamp(heat.colors(128), T_est,min(T_est),max(T_est))
  # plot3d(x = coord[,1], y = coord[,2], z = coord[,3], col = cols,alpha = 0.5,add=TRUE,size=0.6,type="s")
  # 
  # contour3d(brain_templete, level = 50, alpha = 0.1, add = TRUE)
  # 
  # 
  # 
  # 
  # rgl.viewpoint(theta = 0, phi = 5,zoom=0.5)
  # 
  # bgplot3d(image.plot( legend.only=TRUE, zlim= c(min(T_est),max(T_est)), nlevel=128, 
  #                      col=heat.colors(128),horizontal = FALSE))
  # 
  # snapshot3d(paste(image_folder,"/","m_values.png",sep=""))
  # 
  # if(gif==TRUE){
  #   # Save like gif
  #   movie3d(
  #     movie="m_values", 
  #     spin3d( axis = c(0, 0, 1), rpm = 7),
  #     duration = 10, 
  #     dir = image_folder,
  #     type = "gif", 
  #     clean = TRUE
  #   )
  #   
  # }
  
  
  ###################################################################################
  ###End  : To Plot m values ########################################################
  ###################################################################################
  
  
  
  
  ###################################################################################
  ###Begin: To Plot m values ########################################################
  ###################################################################################
  
  
  open3d()
  # Needs to be a bigger window than the default
  par3d(windowRect = c(100*2, 100*2, 1000, 1000))
  
  
  
  
  for(v in 1:(V-1)){
    if(Fiber_Ind[v]==Fiber_Ind[v+1]){
      rgl.linestrips(c(coord[v,1], coord[v+1,1]),
                     c(coord[v,2], coord[v+1,2]),
                     c(coord[v,3], coord[v+1,3]),col="yellow")
    }
    
    if(v%in%selected){
      arrow3d(coord[v,1:3], coord[v,1:3]+ Ref_est[,v]*1.5, type = "lines", col = "Red")
      arrow3d(coord[v,1:3], coord[v,1:3]+ R_est[,v]*1.5, type = "lines", col = "Green")
    }
  }
  
  #cols <- myColorRamp(heat.colors(128), T_est,min(T_est),max(T_est))
  #plot3d(x = coord[,1], y = coord[,2], z = coord[,3], col = cols,alpha = 0.5,add=TRUE,size=0.2,type="s")
  
  contour3d(brain_templete, level = 50, alpha = 0.1, add = TRUE)
  
  
  
  
  rgl.viewpoint(theta = 0, phi = 0,zoom=0.25)
  
  
  
  snapshot3d(paste(image_folder,"/","directional_values.png",sep=""))
  
  if(gif==TRUE){
    # Save like gif
    movie3d(
      movie="directional_values", 
      spin3d( axis = c(0, 0, 1), rpm = 7),
      duration = 10, 
      dir = image_folder,
      type = "gif", 
      clean = TRUE
    )
    
  }
  close3d()
  ###################################################################################
  ###End: To Plot m values ########################################################
  ###################################################################################
  
  
  
}
