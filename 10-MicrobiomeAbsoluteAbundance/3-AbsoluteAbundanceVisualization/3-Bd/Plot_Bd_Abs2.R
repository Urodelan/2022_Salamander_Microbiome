#######################################
### Plot Bd And Bd-Status Densities ###
#######################################

# Better title.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####################
### Begin Script ###
####################

# Read in fungal density predictions.
fungal_densities<-read.csv(file="../RegionITS/FromCluster/RegionITS_NB_Reg_Density_Predictions.csv",stringsAsFactors=F)

# Get just Bd density predictions.
densities<-fungal_densities[fungal_densities$Taxon=="Fungi;Chytridiomycota;Rhizophydiomycetes;Rhizophydiales;Rhizophydiales_fam_Incertae_sedis;Batrachochytrium;dendrobatidis",]

# Read in Bd-status predictions.
status_densities<-read.csv(file="../Region16S/FromCluster/Region16S_NB_Reg_Antifungal_Predictions.csv",stringsAsFactors=F)

# Remove 'other' bacterial taxa from the Bd-status predictions.
status_densities<-status_densities[status_densities$Taxon!="Other",]

# Rename Bd-inhibition status options.
status_densities$Taxon<-ifelse(status_densities$Taxon=="Inhibitory","Bd-Inhibitory Bacterial Taxa",status_densities$Taxon)
status_densities$Taxon<-ifelse(status_densities$Taxon=="NonInhibitory","Non-Bd-Inhibitory Bacterial Taxa",status_densities$Taxon)
status_densities$Taxon<-ifelse(status_densities$Taxon=="Uncertain","Uncertain Bd-Inhibition Status Bacterial Taxa",status_densities$Taxon)

# Combine Bd and Bd-status prediction data frames.
densities<-rbind(densities,status_densities)

# Format date field.
densities$Date<-as.Date(densities$Date)

# Re-format the date field.
densities$Date<-format(densities$Date,"%b %e")

# Remove extra space from the date field.
densities$Date<-gsub("  "," ",densities$Date)

# Format date field as a factor.
densities$Date<-factor(densities$Date,levels=unique(densities$Date))

# Create a numeric date field.
densities$Date_numeric<-as.numeric(densities$Date)

# Define function for getting the most precise taxonomic level.
getMostPreciseTaxonomicLevel<-function(name){
  # Split taxon name by semi-colon.
  x<-strsplit(name,split=";")[[1]]
  # If the taxon is bacterial.
  if(x[1]=="Bacteria"){
    # Store the trailing number part.
    num<-gsub(pattern="_",replacement="",x=x[length(x)])
    # Remove the trailing number part from the name.
    x<-x[-length(x)]
    # If species is included.
    if(length(x)==7){
      # Then get both the binomial name.
      x<-paste(x[6],x[7],num)
    } else {
      # Otherwise, get the most precise taxonomic level available.
      x<-paste(x[length(x)],num)
    }
  } else { # If the taxon is fungal.
    # If species is included.
    if(length(x)==7){
      # Then get both the binomial name.
      x<-paste(x[6],x[7])
    } else {
      # Otherwise, get the most precise taxonomic level available.
      x<-x[length(x)]
    }
  }
  # Return the name.
  return(x)
}

# Create empty field for taxon simple names.
densities$Taxon_simple_name<-NA

# Loop through each record.
for(i in 1:nrow(densities)){
  # Get the taxon's simple name.
  densities$Taxon_simple_name[i]<-getMostPreciseTaxonomicLevel(name=densities$Taxon[i])
}

# Format taxa simple name as a factor to control the ordering in the facet plot.
densities$Taxon_simple_name<-factor(densities$Taxon_simple_name,levels=c("Batrachochytrium dendrobatidis","Bd-Inhibitory Bacterial Taxa","Non-Bd-Inhibitory Bacterial Taxa","Uncertain Bd-Inhibition Status Bacterial Taxa"))

# Create life stage field.
densities$Life_Stage<-paste(densities$Age,ifelse(densities$LM=="Larvae/Neotene",ifelse(densities$Age=="Age-2+","Neotene","Larvae"),"Metamorphosed"))

# Format life stage field as a factor.
densities$Life_Stage<-factor(densities$Life_Stage,levels=c("Age-0 Larvae","Age-0 Metamorphosed","Age-1 Larvae","Age-1 Metamorphosed","Age-2+ Neotene","Age-2+ Metamorphosed"))

# Create data frame storing information on shading by site.
shading<-data.frame(min=seq(from=0.5,to=max(densities$Date_numeric)-0.5,by=1),
                    max=seq(from=1.5,to=max(densities$Date_numeric)+0.5,by=1),
                    col=c("Gibson Lakes","Ponds Lake"))

# Load ggplot2.
library(ggplot2)

# Load package containing lighten function.
library(colorspace)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Load package for creating nice y-axis labels on the log10 scale.
library(scales)

# Set densities plot position dodge value.
densities_position_dodge<-0.75

# Plot predicted densities.
DensitiesPlot<-ggplot(densities,aes(x=Date_numeric,y=Quantile_0.5,color=Life_Stage))+
  facet_wrap(facets="Taxon_simple_name",
             ncol=2,
             nrow=2,
             scales="free",
             strip.position="top",
             dir="h")+
  geom_rect(data=shading,
            aes(xmin=min,xmax=max,ymin=0,ymax=Inf,fill=factor(col)),
            alpha=0.4,
            inherit.aes=F)+
  geom_linerange(aes(ymin=Quantile_0.025,ymax=Quantile_0.975),
                 size=0.5,position=position_dodge(width=densities_position_dodge))+
  geom_linerange(aes(ymin=Quantile_0.25,ymax=Quantile_0.75),
                 size=0.75,position=position_dodge(width=densities_position_dodge))+
  geom_point(size=1,position=position_dodge(width=densities_position_dodge))+
  ggtitle("Densities of Bd-Inhibition Categories on Salamander Skin",
          subtitle="Point: Posterior median \u2013 Thick line: 50% credible interval \u2013 Thin line: 95% credible interval")+
  ylab("Density (arbitrary units)")+
  xlab(NULL)+
  theme_light()+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(),
        panel.grid.minor.y=element_blank(),
        text=element_text(family="Arial"))+
  scale_color_manual(values=
                       c("gold","gold3","green2","green4",lighten("dodgerblue1",0.25),"dodgerblue3"),
                     name="Stage Class")+
  scale_fill_manual(values=c("white","gray53"),name="Site")+
  guides(fill=guide_legend(override.aes=list(color="black",size=0.25)))+
  scale_x_continuous(breaks=unique(densities$Date_numeric),
                     labels=levels(densities$Date),
                     expand=expansion(mult=c(0,0)))+
  scale_y_log10(expand=expansion(mult=c(0.05,0.05)),
                breaks=trans_breaks("log10",function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)))

# Print plot.
print(DensitiesPlot)

# Save densities plot.
ggsave(filename="DensitiesPlot_Bd2.pdf",device=cairo_pdf,plot=DensitiesPlot,width=14,height=10)

# Done!
print("Done!")
