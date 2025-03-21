####### Loading Required Libraries #############################################
library(ggplot2)
library(ggfortify)
library(miceadds)
################################################################################

Base="/Users/zhoulan/Dropbox (Personal)/AkapravaRoy_ZhouLan/Technometrics_Revision/"

####### Defining Paths #########################################################
Output_Path=paste0(Base,"VMF_Regression_Codes/4_1_Method_Comparison/Simulation_Output/")
Fig_Path=paste0(Base,"VMF_Regression_Codes/4_1_Method_Comparison/Figures/")
################################################################################


kappa_all=c(20,30,40)
P_all=2:7
seed_all=1:50
comb_normal=expand.grid(seed_all,kappa_all)
comb_vmf=expand.grid(seed_all,kappa_all,P_all)



RR=data.frame(Method=NULL,alpha=NULL,Metric=NULL,kappa=NULL,Type=NULL)
#### Loading Results of Normal 1 ####
for(index in 1:nrow(comb_normal)){
  
  Rep=comb_normal[index,1]
  kappa=comb_normal[index,2]
  
  load(paste(Output_Path,"VMF_Normal1_","_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
  temp1=data.frame(Method="Gaussian Regression 1",
                   alpha=mean(Validation_Normal1$alpha),
                   Metric="Seperation Angle",Type="Benchmark",kappa=kappa)
  temp2=data.frame(Method="Gaussian Regression 1",
                   alpha=mean(Validation_Normal1.regular$alpha),
                   Metric="Root Mean Squared Error",Type="Benchmark",kappa=kappa)
  RR=rbind(RR,temp1,temp2)
}
###########################################



#### Loading Results of Normal 2 ####
for(index in 1:nrow(comb_normal)){
  
  Rep=comb_normal[index,1]
  kappa=comb_normal[index,2]
  
  load(paste(Output_Path,"VMF_Normal2_","_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
  temp1=data.frame(Method="Gaussian Regression 2",
                   alpha=mean(Validation_Normal2$alpha),
                   Metric="Seperation Angle",Type="Benchmark",kappa=kappa)
  temp2=data.frame(Method="Gaussian Regression 2",
                   alpha=mean(Validation_Normal2.regular$alpha),
                   Metric="Root Mean Squared Error",Type="Benchmark",kappa=kappa)
  RR=rbind(RR,temp1,temp2)
}
############################################

#### Loading Results of NonSpatial VMF ####
for(index in 1:nrow(comb_normal)){
  
  Rep=comb_normal[index,1]
  kappa=comb_normal[index,2]
  
  load(paste(Output_Path,"NonSpatial_VMF_FITTING","_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))
  temp1=data.frame(Method="Non-Spatial VMF",
                   alpha=mean(Validation_lmvmf_p$alpha),
                   Metric="Seperation Angle",Type="Benchmark",kappa=kappa)
  temp2=data.frame(Method="Non-Spatial VMF",
                   alpha=mean(Validation_lmvmf_p.regular$alpha),
                   Metric="Root Mean Squared Error",Type="Benchmark",kappa=kappa)
  RR=rbind(RR,temp1,temp2)
}
############################################

#### Loading Results of VMF ####
for(index in 1:nrow(comb_vmf)){
  
  Rep=comb_vmf[index,1]
  kappa=comb_vmf[index,2]
  P=comb_vmf[index,3]
  
load(paste(Output_Path,"VMF_FITTING_P_",P,"_kappa_",kappa,"_Rep_",Rep,".Rdata",sep=""))

    temp1=data.frame(Method=paste0("VMF: P=",P),
                     alpha=mean(Validation_lmvmf_p$alpha),
                     Metric="Seperation Angle",Type="Proposed",kappa=kappa)
    temp2=data.frame(Method=paste0("VMF: P=",P),
                     alpha=mean(Validation_lmvmf_p.regular$alpha),
                     Metric="Root Mean Squared Error",Type="Proposed",kappa=kappa)
  
  
  
  RR=rbind(RR,temp1,temp2)
}
############################################

pdf(paste0(Fig_Path,"Seperation Angle.pdf"), width = 8, height = 4)
P1=ggplot(RR[RR$Metric=="Seperation Angle",], aes(x=Method, y=alpha, color=Type)) +
  geom_boxplot(outlier.shape = NA)+facet_wrap(~kappa,nrow=1)+
  theme(legend.position="none",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 75, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=12)
  )+ylab("Seperation Angle")
show(P1)
dev.off()


pdf(paste0(Fig_Path,"Root Mean Squared Error.pdf"), width = 8, height = 4)
P2=ggplot(RR[RR$Metric=="Root Mean Squared Error",], aes(x=Method, y=alpha, color=Type)) +
  geom_boxplot(outlier.shape = NA)+ylim(0,1)+facet_wrap(~kappa,nrow=1)+
  theme(legend.position="none",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 75, hjust = 1),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=12)
  )+ylab("Root Mean Squared Error")
show(P2)
dev.off()


