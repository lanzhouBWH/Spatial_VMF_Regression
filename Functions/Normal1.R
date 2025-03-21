Normal1<-function(Data_list,iters=500,burnins=0){
 
  #### Exact List Object to environment#################################################################
  list2env(Data_list, globalenv())
  I=length(Y); V=length(Fiber_Ind); P=ncol(X); K=length(unique(Fiber_Ind)); G=length(unique(group_index))
  ######################################################################################################
  
  
Y_array=array(NA,dim=c(I,V,3))
for(i in 1:I){
  Y_array[i,,]<-Y[[i]]
}



model_string <- nimbleCode({
  
  

  
  for (i in 1:I){
    for(v in 1:V){
      Y_array[i,v,1:3]~dmnorm(mean=MU[i,v,1:3],prec=S[1:3,1:3])
      
      #Y_array[i,v,1]~dnorm(mean=MU[i,v,1],var=S[1])
      #Y_array[i,v,2]~dnorm(mean=MU[i,v,2],var=S[2])
      #Y_array[i,v,3]~dnorm(mean=MU[i,v,3],var=S[3])
      
      MU[i,v,1]<-inprod(X[i,1:P],beta1[group_index[i],v,1:P])
      MU[i,v,2]<-inprod(X[i,1:P],beta2[group_index[i],v,1:P])
      MU[i,v,3]<-inprod(X[i,1:P],beta3[group_index[i],v,1:P])
    }
  }
  
  for (g in 1:G){
    for(v in 1:V){
      beta1[g,v,1:P]~dmnorm(mean=MM[1:P],prec=ZZ[1:P,1:P])
      beta2[g,v,1:P]~dmnorm(mean=MM[1:P],prec=ZZ[1:P,1:P])
      beta3[g,v,1:P]~dmnorm(mean=MM[1:P],prec=ZZ[1:P,1:P])
    }
  }
  
  
  
  S[1:3,1:3]~dwish(RR[1:3,1:3],4)
  
  
})



consts   <- list(I=I,V=V,P=P,G=G,RR=diag(3),MM=rep(0,P),ZZ=diag(P)/100,X=X,group_index=group_index)
data     <- list(Y_array=(Y_array))
inits    <- list(S=diag(3),MU=Y_array)
SS  <- nimbleMCMC(model_string, data = data, inits = inits,
                  constants=consts,
                  samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                  niter = iters, nburnin=burnins , nchains = 1)
return(SS)

}
