####### Loading Required Libraries #############################################
library(miceadds)
################################################################################

#### Please specify the Base folder based on your own computer #################
Base="/Users/zhoulan/Dropbox (Personal)/AkapravaRoy_ZhouLan/Technometrics_Revision/"

####### Defining Paths #########################################################
Data_Path=paste0(Base,"VMF_Regression_Codes/4_1_Method_Comparison/Data/") 
Output_Path=paste0(Base,"VMF_Regression_Codes/4_1_Method_Comparison/Simulation_Output/")
Function_Path=paste0(Base,"VMF_Regression_Codes/Functions/")
################################################################################

####### Source all the functions ###############################################
source.all(Function_Path)
################################################################################

####### Load Synthetic Data ####################################################
load(paste0(Data_Path,"synthetic.Rdata"))
################################################################################

###### Specify Simulation Conditions ###########################################
Rep=1
kappa=40
P=2
set.seed(Rep)
################################################################################


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
fit<-VMF_Reg(Data_List_Training,iters=Iter,burnins=Burn,adaptive=Burn,P=P,Para = TRUE) ### MCMC Model Fitting
################################################################################
################################################################################
################################################################################


################################################################################
######## Step 3: Evaluation Metric Calculation #################################
################################################################################
Validation_lmvmf_p<-Validation.lmvmf(fit,Data_List_Validation,iters=fit$it,burnins=0)
Validation_lmvmf_p.regular<-Validation.lmvmf.regular(fit,Data_List_Validation,iters=fit$it,burnins=0)
save(Validation_lmvmf_p, Validation_lmvmf_p.regular, file=paste(Output_Path,"VMF_FITTING_P_",P,"_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
################################################################################
################################################################################
################################################################################
