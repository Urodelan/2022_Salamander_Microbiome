##################################################
### Fit Model Between Bd-Inhibotry Taxa and Bd ###
##################################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in salamander antifungal table.
sal<-read.csv(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Antifungal_Tables/Salamander_Antifungal_Table.csv",row.names=1,check.names=F)

# Read in salamander fungal taxa table.
bd<-read.csv(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Full_Tables/Region_ITS/Salamander_ITS_Table.csv",row.names=1,check.names=F)

# Get just Bd.
bd<-data.frame(SampleID=row.names(bd),BdReads=bd$`Fungi;Chytridiomycota;Rhizophydiomycetes;Rhizophydiales;Rhizophydiales_fam_Incertae_sedis;Batrachochytrium;dendrobatidis`,TotalBioReadsITS=rowSums(bd),row.names=1:nrow(bd))

# Read in metadata.
meta<-read.csv(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Metadata/SampleMetadata2.csv",stringsAsFactors=F)

# Get the number of salamanders which were stored in the same bucket.
# num_sals_in_buckets<-aggregate(Type~Date+Stratum+Age,data=meta[meta$Type=="Salamander",],FUN=length)
# hist(num_sals_in_buckets$Type)
# summary(num_sals_in_buckets$Type)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 1.000   2.000   3.000   3.044   4.000   5.000
## Mean and median of 3 salamanders stored in a bucket together.

# Add Bd counts to the metadata.
meta<-merge(x=meta,y=bd,all.x=F,all.y=T)

# # Add a field for whether Bd was detected.
# meta$BdDet<-ifelse(meta$BdReads > 0, 1, 0)
# 
# # Subset metadata to just samples in which Bd was detected.
# meta<-meta[meta$BdDet==1,]
# 
# # What percent of samples with Bd detection were Ponds Lake age-2+ metamorphosed salamanders?
# paste0(round(sum(meta$Site=="Ponds Lake" & meta$LM=="Metamorphosed" & meta$Age=="Age-2+")/nrow(meta)*100,1),"%")
# ## 70.7%

# Check how many larval salamander had Bd detected on their skin.
# larvae<-meta[meta$LM=="Larvae/Neotene",]
# sum(larvae$BdReads > 0) # 3
# nrow(larvae) # 139
# Only 3 out of 139 larval salamanders (2.2%) had Bd detected on their skin.

# Subset metadata to just metamorphosed salamanders.
meta<-meta[meta$LM=="Metamorphosed",]
# sum(meta$BdReads > 0) # 38
# nrow(meta) # 66
# 38 out of 66 metamorphosed salamanders (57.6%) had Bd detected on their skin.

# Subset salamander 16S read counts to just samples included in the metadata.
sal<-sal[row.names(sal) %in% meta$SampleID,]

# Subset metadata to just 16S samples.
meta<-meta[meta$SampleID %in% row.names(sal),]

# Re-order salamander 16S data to match the metadata.
sal<-sal[match(meta$SampleID,row.names(sal)),]

# Check that metadata and salamander 16S data are ordered the same.
identical(meta$SampleID,row.names(sal)) # True

# Create a field for total 16S biological 16S reads.
sal$TotalBioReads16S<-rowSums(sal)

# Create fields for the proportion of inhibitory reads and the proportion of Bd.
sal$PropInhib<-sal$Inhibitory/sal$TotalBioReads16S
meta$PropBd<-meta$BdReads/meta$TotalBioReadsITS

# Create a data frame for a beta-binomial regression.
data_bb<-data.frame(SampleID=meta$SampleID,ReadsBd=meta$BdReads,TotalBioReadsITS=meta$TotalBioReadsITS,ReadsInhib=sal$Inhibitory,TotalBioReads16S=sal$TotalBioReads16S,PropBd=meta$PropBd,PropInhib=sal$PropInhib,stringsAsFactors=F)

# Define the logit function.
logit<-function(p) log(p/(1-p))

# Define the inverse logit function.
inv.logit<-function(x) exp(x)/(1+exp(x))

# Plot the logit of proportions against each other.
plot(logit(PropBd)~logit(PropInhib),data=data_bb)

# Create a field for the logit of proportion of inhibitory reads.
data_bb$logit_PropInhib<-logit(data_bb$PropInhib)

# Load rjags.
library(rjags)

# Define model data.
BB_Reg_Data<-list(
  NSamples=nrow(data_bb), # Number of samples.
  BdReadCount=data_bb$ReadsBd, # Bd read counts.
  TotalBioReadsITS=data_bb$TotalBioReadsITS, # Total biological ITS reads.
  InhibReadCount=data_bb$ReadsInhib, # Inhibitory read counts.
  TotalBioReads16S=data_bb$TotalBioReads16S # Total biological 16S reads.
)

# # Initialize the model.
# BB_Reg<-jags.model(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Model_Exp_JAGS_AllMeta.txt",
#                    data=BB_Reg_Data,
#                    n.adapt=20000,
#                    n.chains=3)
# 
# # Burn-in the model.
# update(BB_Reg,n.iter=20000)
# 
# # Sample the model.
# BB_Reg_Out<-coda.samples(model=BB_Reg,
#                          variable.names=c("beta0","betaInhib","phi","PropInhib"),
#                          n.iter=100000,
#                          thin=1)

# Load model output.
load("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Out_AllMeta.RData")

# # Traceplots.
# par(mar=c(1,1,1,1))
# plot(BB_Reg_Out,ask=F)
# par(mar=c(5.1,4.1,4.1,2.1))

# Get convergence summaries.
(gelman_diag_values<-gelman.diag(BB_Reg_Out,multivariate=F)$psrf)

# # Save model output.
# save(BB_Reg_Out,file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Out_AllMeta.RData")

# #### Check traceplots of sample proportions for the following samples.
# ## PropInhib[8]
# ## PropInhib[15]
# ## PropInhib[20]
# 
# # Get numeric target columns.
# target_cols<-which(colnames(BB_Reg_Out[[1]]) %in% c("PropInhib[8]","PropInhib[15]","PropInhib[20]"))
# 
# # Get all chains.
# chain1<-as.data.frame(BB_Reg_Out[[1]])[,target_cols]
# chain2<-as.data.frame(BB_Reg_Out[[2]])[,target_cols]
# chain3<-as.data.frame(BB_Reg_Out[[3]])[,target_cols]
# 
# # Combine the chains into a single data frame.
# chains<-cbind(data.frame(Chain=c(rep(1,nrow(chain1)),rep(2,nrow(chain2)),rep(3,nrow(chain3))),
#                        Draw=c(1:nrow(chain1),1:nrow(chain2),1:nrow(chain3))),
#             rbind(chain1,chain2,chain3))
# 
# # Format the chains field as a factor.
# chains$Chain<-as.factor(chains$Chain)
# 
# # Load ggplot2.
# library(ggplot2)
# 
# # Traceplot of PropInhib[8].
# # Looks fine.
# ggplot()+
#   geom_line(aes(x=Draw,y=`PropInhib[8]`,color=Chain),data=chains)
# 
# # Traceplot of PropInhib[15].
# # Looks fine.
# ggplot()+
#   geom_line(aes(x=Draw,y=`PropInhib[15]`,color=Chain),data=chains)
# 
# # Traceplot of PropInhib[20].
# # Looks alright.
# ggplot()+
#   geom_line(aes(x=Draw,y=`PropInhib[20]`,color=Chain),data=chains)

# Write out model Gelman convergence diagnostics.
write.csv(gelman_diag_values,file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Exp_Convergence_AllMeta.csv",row.names=T)

# # Model summary.
# summary(BB_Reg_Out)

# Combine all MCMC estimates.
BB_Reg_MCMC<-as.data.frame(do.call("rbind",BB_Reg_Out))

# What is the probability that there is a negative relationship between the proportion
# of antifungal reads and the proportion of Bd reads?
paste0(round(mean(BB_Reg_MCMC$betaInhib < 0)*100,1),"%") # 100%

# Get summary statistics of MCMC samples.
(BB_Reg_MCMC_quantiles<-as.data.frame(t(apply(BB_Reg_MCMC,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975)))))

# Write out summary statistics of MCMC samples.
write.csv(BB_Reg_MCMC_quantiles,file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Exp_MCMCQuantiles_AllMeta.csv",row.names=T)

# Create a vector of predictor values.
pred_val<-seq(from=min(data_bb$logit_PropInhib),to=max(data_bb$logit_PropInhib),length.out=1e3)

# This method uses too much memory.
# # Generate posterior predictions.
# BB_predictions<-BB_Reg_MCMC$beta0+as.matrix(BB_Reg_MCMC$betaInhib) %*% t(pred_val)
# 
# # Get summaries of the posterior predictions.
# BB_prediction_quantiles<-as.data.frame(t(apply(X=BB_predictions,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
# 
# # Rename fields in the posterior prediction summary data frame.
# colnames(BB_prediction_quantiles)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=colnames(BB_prediction_quantiles)))/100)

# Create empty storage data frame for prediction quantiles.
BB_prediction_quantiles_df<-data.frame(NULL)

# Loop through each predictor values.
for(i in 1:length(pred_val)){
  
  # Generate predictions from the MCMC samples.
  BB_prediction<-BB_Reg_MCMC$beta0+BB_Reg_MCMC$betaInhib*pred_val[i]
  
  # Get prediction quantiles.
  BB_prediction_quantiles<-quantile(x=BB_prediction,probs=c(0.025,0.25,0.5,0.75,0.975))
  
  # Rename elements of the prediction quantile vector.
  names(BB_prediction_quantiles)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=names(BB_prediction_quantiles)))/100)
  
  # Create a data frame containing prediction quantiles and the predictor value.
  BB_prediction_quantiles<-cbind(data.frame(logit_PropInhib=pred_val[i]),
                                 as.data.frame(t(BB_prediction_quantiles)))
  
  # Add the prediction quantiles and predictor value to the storage data frame.
  BB_prediction_quantiles_df<-rbind(BB_prediction_quantiles_df,BB_prediction_quantiles)

}

# Load ggplot2.
library(ggplot2)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Create plot of posterior predictions.
BB_Reg_logit_plot<-ggplot()+
  geom_point(aes(x=logit_PropInhib,y=logit(PropBd)),data=data_bb)+
  geom_ribbon(aes(x=logit_PropInhib,ymin=Quantile_0.025,ymax=Quantile_0.975),data=BB_prediction_quantiles_df,fill=NA,color="black",linetype="dashed")+
  geom_line(aes(x=logit_PropInhib,y=Quantile_0.5),data=BB_prediction_quantiles_df,color="black",linetype="solid")+
  theme_light()+
  xlab(expression("logit(Proportional Abundance of"~italic("Bd")*"-Inhibitory Bacterial Taxa)"))+
  ylab(expression("logit(Proportional Abundance of"~italic("Bd")*")"))+
  ggtitle("Bayesian Beta-Binomial Regression for Metamorphosed Salamanders",subtitle="Point: Observation \u2013 Solid line: Posterior median \u2013 Dashed line: 95% credible interval")+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        text=element_text(family="Arial"))+
  scale_y_continuous(expand=expansion(mult=0.025))+
  scale_x_continuous(expand=expansion(mult=0.025))

# Display plot.
print(BB_Reg_logit_plot)

# Save plot.
ggsave(filename="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Exp_logit_AllMeta2.pdf",device=cairo_pdf,plot=BB_Reg_logit_plot,width=7,height=5)

# Now generate predictions on the proportions scale.

# Create a vector of predictor values.
pred_val<-seq(from=min(data_bb$PropInhib),to=max(data_bb$PropInhib),length.out=1e3)

# Create empty storage data frame for prediction quantiles.
BB_prediction_quantiles_df<-data.frame(NULL)

# Loop through each predictor values.
for(i in 1:length(pred_val)){
  
  # Generate predictions from the MCMC samples.
  BB_prediction<-inv.logit(BB_Reg_MCMC$beta0+BB_Reg_MCMC$betaInhib*logit(pred_val[i]))
  
  # Get prediction quantiles.
  BB_prediction_quantiles<-quantile(x=BB_prediction,probs=c(0.025,0.25,0.5,0.75,0.975))
  
  # Rename elements of the prediction quantile vector.
  names(BB_prediction_quantiles)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=names(BB_prediction_quantiles)))/100)
  
  # Create a data frame containing prediction quantiles and the predictor value.
  BB_prediction_quantiles<-cbind(data.frame(PropInhib=pred_val[i]),
                                 as.data.frame(t(BB_prediction_quantiles)))
  
  # Add the prediction quantiles and predictor value to the storage data frame.
  BB_prediction_quantiles_df<-rbind(BB_prediction_quantiles_df,BB_prediction_quantiles)
  
}

# Create plot of posterior predictions.
BB_Reg_standard_plot<-ggplot()+
  geom_point(aes(x=PropInhib,y=PropBd),data=data_bb)+
  geom_ribbon(aes(x=PropInhib,ymin=Quantile_0.025,ymax=Quantile_0.975),data=BB_prediction_quantiles_df,fill=NA,color="black",linetype="dashed")+
  geom_line(aes(x=PropInhib,y=Quantile_0.5),data=BB_prediction_quantiles_df,color="black",linetype="solid")+
  theme_light()+
  xlab(expression("Proportional Abundance of"~italic("Bd")*"-Inhibitory Bacterial Taxa"))+
  ylab(expression("Proportional Abundance of"~italic("Bd")))+
  ggtitle("Bayesian Beta-Binomial Regression for Metamorphosed Salamanders",subtitle="Point: Observation \u2013 Solid line: Posterior median \u2013 Dashed line: 95% credible interval")+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        text=element_text(family="Arial"))+
  scale_y_continuous(expand=expansion(mult=c(0.01,0.01)))+
  scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(NA,0.01)))

# Display plot.
print(BB_Reg_standard_plot)

# Save plot.
ggsave(filename="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Exp_standard_AllMeta2.pdf",device=cairo_pdf,plot=BB_Reg_standard_plot,width=7,height=5)

# Remove title and subtitle from plots.
BB_Reg_logit_plot<-BB_Reg_logit_plot+ggtitle(NULL,subtitle=NULL)
BB_Reg_standard_plot<-BB_Reg_standard_plot+ggtitle(NULL,subtitle=NULL)

# Load grid.
library(grid)

# Create grob object.
grob1<-grobTree(textGrob("logit scale",x=0.98,y=0.98,hjust=1,vjust=1,
                          gp=gpar(col="black",fontsize=14,fontface="italic")))

# Add grob object to plot.
bb_logit_text<-BB_Reg_logit_plot+annotation_custom(grob1)

# Create grob object.
grob2<-grobTree(textGrob("Proportions",x=0.98,y=0.98,hjust=1,vjust=1,
                         gp=gpar(col="black",fontsize=14,fontface="italic")))

# Add grob object to plot.
bb_standard_text<-BB_Reg_standard_plot+annotation_custom(grob2)

# Load cowplot.
library(cowplot)

# Combine plots.
combined_plots<-plot_grid(bb_logit_text,bb_standard_text,ncol=2)

# Make a graphical supertitle object.
super_title<-ggdraw()+
  draw_label("Bayesian Beta-Binomial Regression for Metamorphosed Salamanders",fontface="bold",hjust=0.5,x=0.5,y=0.67)+
  draw_label("Point: Observation \u2013 Solid line: Posterior median \u2013 Dashed line: 95% credible interval",x=0.5,y=0.25)+
  theme(text=element_text(family="Arial"))

# Add supertitle to the plot panel.
combined_plots_full<-plot_grid(super_title,combined_plots,ncol=1,rel_heights=c(0.125,1))

# Save plot.
ggsave(filename="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/Full_Taxa_Tables/Bd_Models/Proportions_AllMeta/BB_Reg_Exp_both_AllMeta2.pdf",device=cairo_pdf,plot=combined_plots_full,width=14,height=6)

# Done!
print("Done!")
