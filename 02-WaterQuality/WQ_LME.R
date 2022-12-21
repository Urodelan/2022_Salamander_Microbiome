###############################################
### Water Quality Data Mixed Effects Models ###
###############################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory and load packages.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in water quality data.
wq<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Field/Field_data/Data/Data/Original_data/CSV/Water_quality.csv",stringsAsFactors=F)

# Format date.
wq$Date<-as.Date(wq$Day,format='%m/%d/%Y')

# Convert date to number of weeks since June 9th, 2018.
wq$Week<-as.vector(wq$Date-as.Date("2018-06-09"))/7

# Ensure that Site is a factor.
wq$Site<-as.factor(wq$Site)

# Make strata for Ponds Lake numbers 5 through 9.
wq$Stratum<-ifelse(wq$Site=="Ponds Lake",wq$Stratum+4,wq$Stratum)

# Load package for fitting linear-mixed effects models with p-values. 
library(lmerTest)

# Get the package version.
# packageVersion("lmerTest") # 3.1.3

# lmerTest citation:
## Kuznetsova A, Brockhoff PB, Christensen RHB (2017). “lmerTest Package: Tests in Linear
## Mixed Effects Models.” _Journal of Statistical Software_, *82*(13), 1-26.
## doi:10.18637/jss.v082.i13 <https://doi.org/10.18637/jss.v082.i13>.

# Regression coefficients with p < 0.05 are listed (rounded values are in parentheses).

# Fit model for temperature.
temp<-lmer(Temperature_C~(1|Stratum)+Site*Week,data=wq)
summary(temp)
## Week: beta = -0.71285 (-0.713), p-value = 2.56e-12 (< 0.001)

# Fit model for pH.
pH<-lmer(pH~(1|Stratum)+Site*Week,data=wq)
summary(pH)
## Site: beta = -1.04773 (-1.048), p-value = 5.37e-05 (< 0.001)
## Week: beta = 0.05310 (0.053), p-value = 0.000472 (< 0.001)

# Fit model for conductivity.
cond<-lmer(Conductivity_uS~(1|Stratum)+Site*Week,data=wq)
summary(cond)
## Site: beta = -74.8220 (-74.822), p-value = 3.08e-13 (< 0.001)
## Week: beta = 2.1972 (2.197), p-value = 0.00013 (< 0.001)
## Interaction: beta = -2.1205 (-2.120), p-value = 0.00980 (0.010)

# Fit model for dissolved oxygen (ppm).
DO_ppm<-lmer(Dissolved_oxygen_ppm~(1|Stratum)+Site*Week,data=wq)
summary(DO_ppm)
## Site: beta = -2.23306 (-2.233), p-value = 0.000819 (0.001)
## Week: beta = -0.12587 (-0.126), p-value = 0.000185 (< 0.001)
## Interaction: beta = 0.15066 (0.151), p-value = 0.002023 (0.002)

# Fit model for dissolved oxygen (percent).
DO_pct<-lmer(Dissolved_oxygen_pct~(1|Stratum)+Site*Week,data=wq)
summary(DO_pct)
## Site: beta = -23.8841 (-23.884), p-value = 0.00639 (0.006)
## Week: beta = -2.8353 (-2.835), p-value = 1.39e-09 (< 0.001)
## Interaction: beta = 1.9007 (1.901), p-value = 0.00159 (0.002)

# Format model output.

# Temperature.
temp_coef<-as.data.frame(summary(temp)$coefficients)
temp_sig<-temp_coef$`Pr(>|t|)` < 0.05
temp_coef$`Pr(>|t|)`<-ifelse(temp_coef$`Pr(>|t|)` < 0.001,"< 0.001",round(temp_coef$`Pr(>|t|)`,3))
temp_coef$`Pr(>|t|)`<-ifelse(temp_sig,paste0(temp_coef$`Pr(>|t|)`,"*"),temp_coef$`Pr(>|t|)`)
temp_coef[,1:4]<-round(temp_coef[,1:4],3)
temp_coef$Model_Parameter<-row.names(temp_coef)
row.names(temp_coef)<-1:nrow(temp_coef)
temp_coef$Water_Quality_Parameter<-"Temperature"

# pH.
pH_coef<-as.data.frame(summary(pH)$coefficients)
pH_sig<-pH_coef$`Pr(>|t|)` < 0.05
pH_coef$`Pr(>|t|)`<-ifelse(pH_coef$`Pr(>|t|)` < 0.001,"< 0.001",round(pH_coef$`Pr(>|t|)`,3))
pH_coef$`Pr(>|t|)`<-ifelse(pH_sig,paste0(pH_coef$`Pr(>|t|)`,"*"),pH_coef$`Pr(>|t|)`)
pH_coef[,1:4]<-round(pH_coef[,1:4],3)
pH_coef$Model_Parameter<-row.names(pH_coef)
row.names(pH_coef)<-1:nrow(pH_coef)
pH_coef$Water_Quality_Parameter<-"pH"

# Conductivity.
cond_coef<-as.data.frame(summary(cond)$coefficients)
cond_sig<-cond_coef$`Pr(>|t|)` < 0.05
cond_coef$`Pr(>|t|)`<-ifelse(cond_coef$`Pr(>|t|)` < 0.001,"< 0.001",round(cond_coef$`Pr(>|t|)`,3))
cond_coef$`Pr(>|t|)`<-ifelse(cond_sig,paste0(cond_coef$`Pr(>|t|)`,"*"),cond_coef$`Pr(>|t|)`)
cond_coef[,1:4]<-round(cond_coef[,1:4],3)
cond_coef$Model_Parameter<-row.names(cond_coef)
row.names(cond_coef)<-1:nrow(cond_coef)
cond_coef$Water_Quality_Parameter<-"Conductivity"

# Dissolved oxygen (ppm).
DO_ppm_coef<-as.data.frame(summary(DO_ppm)$coefficients)
DO_ppm_sig<-DO_ppm_coef$`Pr(>|t|)` < 0.05
DO_ppm_coef$`Pr(>|t|)`<-ifelse(DO_ppm_coef$`Pr(>|t|)` < 0.001,"< 0.001",round(DO_ppm_coef$`Pr(>|t|)`,3))
DO_ppm_coef$`Pr(>|t|)`<-ifelse(DO_ppm_sig,paste0(DO_ppm_coef$`Pr(>|t|)`,"*"),DO_ppm_coef$`Pr(>|t|)`)
DO_ppm_coef[,1:4]<-round(DO_ppm_coef[,1:4],3)
DO_ppm_coef$Model_Parameter<-row.names(DO_ppm_coef)
row.names(DO_ppm_coef)<-1:nrow(DO_ppm_coef)
DO_ppm_coef$Water_Quality_Parameter<-"Dissolved Oxygen (ppm)"

# Dissolved oxygen (%).
DO_pct_coef<-as.data.frame(summary(DO_pct)$coefficients)
DO_pct_sig<-DO_pct_coef$`Pr(>|t|)` < 0.05
DO_pct_coef$`Pr(>|t|)`<-ifelse(DO_pct_coef$`Pr(>|t|)` < 0.001,"< 0.001",round(DO_pct_coef$`Pr(>|t|)`,3))
DO_pct_coef$`Pr(>|t|)`<-ifelse(DO_pct_sig,paste0(DO_pct_coef$`Pr(>|t|)`,"*"),DO_pct_coef$`Pr(>|t|)`)
DO_pct_coef[,1:4]<-round(DO_pct_coef[,1:4],3)
DO_pct_coef$Model_Parameter<-row.names(DO_pct_coef)
row.names(DO_pct_coef)<-1:nrow(DO_pct_coef)
DO_pct_coef$Water_Quality_Parameter<-"Dissolved Oxygen (%)"

# Combine model output.
wq_coef<-rbind(temp_coef,pH_coef,cond_coef,DO_ppm_coef,DO_pct_coef)

# Re-order fields.
wq_coef<-wq_coef[,c(7,6,1:5)]

# Rename fields.
colnames(wq_coef)<-c("Water Quality Parameter","Model Parameter","Regression Coefficient","Standard Error","Degrees of Freedom","t-statistic","p-value")

# Rename model parameters.
## Intercept term.
wq_coef$`Model Parameter`<-ifelse(wq_coef$`Model Parameter`=="(Intercept)","Intercept",wq_coef$`Model Parameter`)
## Site term.
wq_coef$`Model Parameter`<-ifelse(wq_coef$`Model Parameter`=="SitePonds Lake","Site",wq_coef$`Model Parameter`)
## Interaction term.
wq_coef$`Model Parameter`<-ifelse(wq_coef$`Model Parameter`=="SitePonds Lake:Week","Interaction",wq_coef$`Model Parameter`)

# Write out model output.
write.csv(x=wq_coef,file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Results/Field/WQ_Models_Output.csv",row.names=F)

# Done!
