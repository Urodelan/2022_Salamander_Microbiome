########################################
### Plot 16S Proportional Abundances ###
########################################

# The only difference between this script and the first one is that
# this script saves plots as pdf files.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####################
### Begin Script ###
####################

# Read in proportion predictions.
props<-read.csv(file="../FromCluster/Predictions/props_summary_16S.csv",stringsAsFactors=F)

# Format date field.
props$Date<-as.Date(props$Date)

# Re-format the date field.
props$Date<-format(props$Date,"%b %e")

# Remove extra space from the date field.
props$Date<-gsub("  "," ",props$Date)

# Format date field as a factor.
props$Date<-factor(props$Date,levels=unique(props$Date))

# Create a numeric date field.
props$Date_numeric<-as.numeric(props$Date)

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
props$Taxon_simple_name<-NA

# Loop through each record.
for(i in 1:nrow(props)){
  # Get the taxon's simple name.
  props$Taxon_simple_name[i]<-getMostPreciseTaxonomicLevel(name=props$Taxon[i])
}

# Format taxa simple name as a factor to preserve the ordering in a facet plot.
props$Taxon_simple_name<-factor(props$Taxon_simple_name,levels=unique(props$Taxon_simple_name))

# Get just salamander taxa.
props_sal<-props[props$Type=="Salamander",]

# Get just environmental taxa.
props_envi<-props[props$Type!="Salamander",]

# Create life stage field.
props_sal$Life_Stage<-paste(props_sal$Age,ifelse(props_sal$LM=="Larvae/Neotene",ifelse(props_sal$Age=="Age-2+","Neotene","Larvae"),"Metamorphosed"))

# Format life stage field as a factor.
props_sal$Life_Stage<-factor(props_sal$Life_Stage,levels=c("Age-0 Larvae","Age-0 Metamorphosed","Age-1 Larvae","Age-1 Metamorphosed","Age-2+ Neotene","Age-2+ Metamorphosed"))

# Get the minimum of environmental 0.025 quantiles for each taxon and date.
props_envi_limits_min<-aggregate(formula=Quantile_0.025~Taxon_simple_name+Date_numeric,
                                 data=props_envi,FUN=min)

# Get the maximum of environmental 0.975 quantiles for each taxon and date.
props_envi_limits_max<-aggregate(formula=Quantile_0.975~Taxon_simple_name+Date_numeric,
                                 data=props_envi,FUN=max)

# Combine minimum and maximum environmental proportional abundances into a single dataframe.
props_envi_limits<-merge(props_envi_limits_min,props_envi_limits_max)

# Create an x-min and x-max field for plotting purposes.
props_envi_limits$x_min<-props_envi_limits$Date_numeric-0.5
props_envi_limits$x_max<-props_envi_limits$Date_numeric+0.5

# Create data frame storing information on shading by site.
shading<-data.frame(min=seq(from=0.5,to=max(props_sal$Date_numeric)-0.5,by=1),
                    max=seq(from=1.5,to=max(props_sal$Date_numeric)+0.5,by=1),
                    col=c("Gibson Lakes","Ponds Lake"))

# Load ggplot2.
library(ggplot2)

# Load package containing lighten function.
library(colorspace)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Get a vector of unique taxa.
unique_taxa<-unique(props_sal$Taxon_simple_name)

# Create a dataframe defining the taxa subsets to include in the same plot panel.
unique_taxa<-data.frame(Taxon=unique_taxa,Subset=rep(1:5,each=21)[1:length(unique_taxa)],
                        stringsAsFactors=F)

# Loop through each plot panel taxa subset.
for(i in 1:max(unique_taxa$Subset)){
  
  # Get taxa names to plot.
  taxa_to_plot<-unique_taxa$Taxon[unique_taxa$Subset==i]
  
  # Subset salamander proportions to the taxa to plot.
  props_sal_taxa_to_plot<-props_sal[props_sal$Taxon_simple_name %in% taxa_to_plot,]
  
  # Subset environmental proportions to the taxa to plot.
  props_envi_limits_to_plot<-props_envi_limits[props_envi_limits$Taxon_simple_name %in% taxa_to_plot,]
  
  # Set proportions plot position dodge value.
  props_position_dodge<-0.75
  
  # Plot predicted proportions.
  PropsPlot<-ggplot(props_sal_taxa_to_plot,aes(x=Date_numeric,y=Quantile_0.5,color=Life_Stage))+
    facet_wrap(facets="Taxon_simple_name",
               ncol=3,
               nrow=7,
               scales="free",
               strip.position="top",
               dir="h")+
    geom_rect(data=shading,
              aes(xmin=min,xmax=max,ymin=0,ymax=Inf,fill=factor(col)),
              alpha=0.4,
              inherit.aes=F)+
    geom_rect(data=props_envi_limits_to_plot,
              aes(xmin=x_min,xmax=x_max,ymin=Quantile_0.025,ymax=Quantile_0.975),
              fill="red",
              alpha=0.125,
              inherit.aes=F)+
    geom_linerange(aes(ymin=Quantile_0.025,ymax=Quantile_0.975),
                   size=0.5,position=position_dodge(width=props_position_dodge))+
    geom_linerange(aes(ymin=Quantile_0.25,ymax=Quantile_0.75),
                   size=0.75,position=position_dodge(width=props_position_dodge))+
    geom_point(size=1,position=position_dodge(width=props_position_dodge))+
    ggtitle(paste0("Bacterial Proportional Abundances on Salamander Skin (Plot ",i,
                   " of ",max(unique_taxa$Subset),")"),
            subtitle="Point: Posterior median \u2013 Thick line: 50% credible interval \u2013 Thin line: 95% credible interval\nShaded red range: Environmental proportional abundance")+
    ylab("Proportional Abundance")+
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
    scale_x_continuous(breaks=unique(props_sal$Date_numeric),
                       labels=levels(props_sal$Date),
                       expand=expansion(mult=c(0,0)))+
    scale_y_sqrt(limits=c(0,NA),expand=expansion(mult=c(0,0.05)))
  
  # Save proportions plot.
  ggsave(filename=paste0("PropsPlot_16S_",i,".pdf"),device=cairo_pdf,plot=PropsPlot,width=21,height=21)
  
}

# Done!
print("Done!")
