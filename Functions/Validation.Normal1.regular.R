Validation.Normal1.regular<-function(fit,Data_list,iters,burnins){
  
  #### Exact List Object to environment#################################################################
  list2env(Data_list, globalenv())
  I=length(Y); V=length(Fiber_Ind); p=ncol(X); K=length(unique(Fiber_Ind)); G=length(unique(group_index))
  ######################################################################################################
  
  
  
  pred_Matrix=array(NA, dim=c(3,I,V,iters-burnins))
  
  fit=fit[,10:ncol(fit)]
  data=array(fit,dim=c(iters,G,V,p,3))
  
  
  index=0
  for(it in (burnins+1):iters ){
    index=index+1
    
    for(i in 1:I){
      for(v in 1:V){
        
        # S=diag(3)
        # S[1,1]<-fit[it,"S[1, 1]"]
        # S[1,2]<-fit[it,"S[1, 2]"]
        # S[1,3]<-fit[it,"S[1, 3]"]
        # S[2,1]<-fit[it,"S[2, 1]"]
        # S[2,2]<-fit[it,"S[2, 2]"]
        # S[2,3]<-fit[it,"S[2, 3]"]
        # S[3,1]<-fit[it,"S[3, 1]"]
        # S[3,2]<-fit[it,"S[3, 2]"]
        # S[3,3]<-fit[it,"S[3, 3]"]
        # 
        
        
        #beta1=sapply(1:p, function(pp) fit[it,paste("beta1[", group_index[i],", ",v,", ", pp,"]",sep="" )] )
        #beta2=sapply(1:p, function(pp) fit[it,paste("beta2[", group_index[i],", ",v,", ", pp,"]",sep="" )] )
        #beta3=sapply(1:p, function(pp) fit[it,paste("beta3[", group_index[i],", ",v,", ", pp,"]",sep="" )] )
        
        beta1=data[it,group_index[i],v,,1]
        beta2=data[it,group_index[i],v,,2]
        beta3=data[it,group_index[i],v,,3]
        
        
        mu=c(X[i,]%*%beta1,X[i,]%*%beta2,X[i,]%*%beta3)
        
        pred_Matrix[,i,v,index]<-mu
        
        
      }
      
    }
    
  }
  
  
  PP_Matrix=array(NA, dim=c(3,I,V))
  
  for(i in 1:I){
    PP_Matrix[,i,]<-sapply(1:V, function(v) rowMeans(pred_Matrix[,i,v,])/sum(rowMeans(pred_Matrix[,i,v,])^2)) 
  }
  
  
  alpha=sapply(1:I, function(i) sapply(1:V, function(v) mean((PP_Matrix[,i,v]-Y[[i]][v,])^2) ))
  # S=diag(3)
  # S[1,1]<-mean(fit[(burnins+1):iters,"S[1, 1]"])
  # S[2,2]<-mean(fit[(burnins+1):iters,"S[2, 2]"])
  # S[3,3]<-mean(fit[(burnins+1):iters,"S[3, 3]"])
  # Log_Lik=sapply(1:I, function(i) sapply(1:V, function(v) dmvnorm(Y[[i]][v,],PP_Matrix[,i,v],solve(S),logged = TRUE)))   
  # 
  RR_Normal1=list(alpha=alpha)
  
  return(RR_Normal1)
  
}