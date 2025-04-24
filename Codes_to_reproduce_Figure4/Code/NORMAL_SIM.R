####### Loading Required Libraries #############################################
library(miceadds)
################################################################################


####### Defining Paths #########################################################
Data_Path=paste0("Spatial_VMF_Regression/Codes_to_reproduce_Figure4/Data/") 
Output_Path=paste0("Spatial_VMF_Regression/Codes_to_reproduce_Figure4/Simulation_Output/")
Function_Path=paste0("Spatial_VMF_Regression/Functions/")
################################################################################

####### Source all the functions ###############################################
source.all(Function_Path)
################################################################################

####### Load Synthetic Data ####################################################
load(paste0(Data_Path,"synthetic.Rdata"))
################################################################################

for(Rep in 1:50){            ### The specification for 50 Replications ###
  set.seed(Rep)              ### Set Seed
  for(kappa in c(20,30,40)){ ### The specification for Conerntration Parameter kappa=20,30,40 ###
      cat("Gaussian Regression Running Rep =", Rep, ", kappa =", kappa,  "\n")
      ################################################################################
      ######## Step 1: Data Simulation ###############################################
      ################################################################################
      
      ### Getting Coefficients and Training Data #####################
      Training<-Training[seq(1,80,4),]
      group_index_cat=(Training$Group)
      group_index=rep(NA,length(group_index_cat))
      group_index[group_index_cat=="CN"]=1
      group_index[group_index_cat=="EMCI"]=2
      group_index[group_index_cat=="LMCI"]=3
      group_index[group_index_cat=="AD"]=4
      X=model.matrix(~as.factor(Sex)+ Age+ MMSCORE+ APGEN,Training)
      Data_List_Training=Synthetic_Data(kappa,X,group_index,COEF,Fiber_Ind,5)
      ##############################################
      
      
      ### Getting Validation Data #####################
      Validation<-Validation[seq(1,80,4),]
      group_index_cat=(Validation$Group)
      group_index=rep(NA,length(group_index_cat))
      group_index[group_index_cat=="CN"]=1
      group_index[group_index_cat=="EMCI"]=2
      group_index[group_index_cat=="LMCI"]=3
      group_index[group_index_cat=="AD"]=4
      X=model.matrix(~as.factor(Sex)+ Age+ MMSCORE+ APGEN,Validation)
      Data_List_Validation=Synthetic_Data(kappa,X,group_index,COEF,Fiber_Ind,5)
      ##############################################
      
      ################################################################################
      ################################################################################
      ################################################################################
      
      
      ################################################################################
      ######## Step 2: MCMC Fitting ##################################################
      ################################################################################
      Iter=5000
      Burn=2000
      fit_Normal1<-Normal1(Data_List_Training,iters=Iter,burnins=Burn) ### MCMC algorithm for Gaussian Regression 1
      fit_Normal2<-Normal2(Data_List_Training,iters=Iter,burnins=Burn) ### MCMC algorithm for Gaussian Regression 2
      ################################################################################
      ################################################################################
      ################################################################################
      
      
      ################################################################################
      ######## Step 3: Evaluation Metric Calculation #################################
      ################################################################################
      Validation_Normal1<-Validation.Normal1(fit_Normal1,Data_List_Validation,iters=5000,burnins=0)
      Validation_Normal1.regular<-Validation.Normal2.regular(fit_Normal1,Data_List_Validation,iters=5000,burnins=0)
      save(Validation_Normal1, Validation_Normal1.regular,file=paste(Output_Path,"VMF_Normal1_","_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
      
      Validation_Normal2<-Validation.Normal2(fit_Normal2,Data_List_Validation,iters=5000,burnins=0)
      Validation_Normal2.regular<-Validation.Normal2.regular(fit_Normal2,Data_List_Validation,iters=5000,burnins=0)
      save(Validation_Normal2, Validation_Normal2.regular,file=paste(Output_Path,"VMF_Normal2_","_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
      ################################################################################
      ################################################################################
      ################################################################################
  }
}
