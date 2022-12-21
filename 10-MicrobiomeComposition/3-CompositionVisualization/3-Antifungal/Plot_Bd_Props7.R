#####################################################
### Plot Bd And Bd-Status Proportional Abundances ###
#####################################################

# This script attempts to make some facet titles italicized
# and some plain font face, but it doesn't work. Facet titles are either
# all plain or all italic. The following warning message is thrown:

# Warning message:
#   Vectorized input to `element_text()` is not officially supported.
# Results may be unexpected or may change in future versions of ggplot2. 

# Bd is italicized in the plot title.

# This script plots the output from the updated antifungal DM regression.

# The background grey and environmental box are also adjusted to be easier to read.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####################
### Begin Script ###
####################

# Read in fungal proportion predictions.
fungal_props<-read.csv(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Results/Composition/RegionITS/FromCluster/Predictions/props_summary_ITS.csv",stringsAsFactors=F)

# Get just Bd proportion predictions.
props<-fungal_props[fungal_props$Taxon=="Fungi;Chytridiomycota;Rhizophydiomycetes;Rhizophydiales;Rhizophydiales_fam_Incertae_sedis;Batrachochytrium;dendrobatidis",]

# Read in Bd-status predictions.
status_props<-read.csv(file="/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Results/Composition/Antifungal/FromCluster/Predictions/props_summary_Antifungal.csv",stringsAsFactors=F)

# Rename Bd-inhibition status options.
status_props$Taxon<-ifelse(status_props$Taxon=="Inhibitory","Bd-Inhibitory Bacterial Taxa",status_props$Taxon)
status_props$Taxon<-ifelse(status_props$Taxon=="NonInhibitory","Non-Bd-Inhibitory Bacterial Taxa",status_props$Taxon)
status_props$Taxon<-ifelse(status_props$Taxon=="Uncertain","Uncertain Bd-Inhibition Status Bacterial Taxa",status_props$Taxon)

# Combine Bd and Bd-status prediction data frames.
props<-rbind(props,status_props)

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

# Format taxa simple name as a factor to control the ordering in the facet plot.
props$Taxon_simple_name<-factor(props$Taxon_simple_name,levels=c("Batrachochytrium dendrobatidis","Bd-Inhibitory Bacterial Taxa","Non-Bd-Inhibitory Bacterial Taxa","Uncertain Bd-Inhibition Status Bacterial Taxa"))

# Get just salamander taxa.
props_sal<-props[props$Type=="Salamander",]

# Get just environmental taxa.
props_envi<-props[props$Type!="Salamander",]

# Create life stage field.
props_sal$Life_Stage<-paste(props_sal$Age,ifelse(props_sal$LM=="Larvae/Neotene",ifelse(props_sal$Age=="Age-2+","Neotene","Larvae"),"Metamorphosed"))

# Format life stage field as a factor.
props_sal$Life_Stage<-factor(props_sal$Life_Stage,levels=c("Age-0 Larvae","Age-0 Metamorphosed","Age-1 Larvae","Age-1 Metamorphosed","Age-2+ Neotene","Age-2+ Metamorphosed"))

# Get the minimum of environmental 0.025 quantiles for each taxon and date.
props_envi_limits_min<-aggregate(Quantile_0.025~Taxon_simple_name+Date_numeric,
                                 data=props_envi,FUN=min)

# Get the maximum of environmental 0.975 quantiles for each taxon and date.
props_envi_limits_max<-aggregate(Quantile_0.975~Taxon_simple_name+Date_numeric,
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

# Copy environmental proportions data frame.
props_envi_limits_to_plot<-props_envi_limits

# Set dash spacing.
spacing<-0.25

# Set dash slope.
slope<-0.2

# Create an empty storage data frame for the slashes.
slashes<-data.frame(NULL)

# Loop through each environmental proportional abundance box.
for(i in 1:nrow(props_envi_limits_to_plot)){
  
  # Get the box boundaries.
  record<-props_envi_limits_to_plot[i,]
  
  # Square-root the min and max y coordinates for backtransforming later.
  record$Quantile_0.025<-sqrt(record$Quantile_0.025)
  record$Quantile_0.975<-sqrt(record$Quantile_0.975)
  
  # Get initial dash coordinates.
  diag<-data.frame(x_min=seq(from=record$x_min-(record$Quantile_0.975-record$Quantile_0.025)/slope,
                             to=record$x_max,by=spacing),
                   y_min=record$Quantile_0.025,
                   y_max=record$Quantile_0.975)
  
  # Derive maximum y and y coordinates.
  diag$x_max<-(diag$y_max-diag$y_min)/slope+diag$x_min
  diag$y_max<-(diag$x_max-diag$x_min)*slope+diag$y_min
  
  # Loop through each dash.
  for(j in 1:nrow(diag)){
    
    # Get the dash coordinates.
    diag_record<-diag[j,]
    
    # If the dash begins before the min x value.
    if(diag_record$x_min < record$x_min){
      
      # Interpolate to create a new starting x and y value.
      diag_record$y_min<-diag_record$y_max-(diag_record$x_max-record$x_min)*slope
      diag_record$x_min<-record$x_min
      
    }
    
    # If the dash ends after the max x value.
    if(diag_record$x_max > record$x_max){
      
      # Interpolate to create a new ending x and y value.
      diag_record$y_max<-diag_record$y_min+(record$x_max-diag_record$x_min)*slope
      diag_record$x_max<-record$x_max
      
    }
    
    # Update the slash coordinates.
    diag[j,]<-diag_record
    
  }
  
  # Add taxon name to the slashes.
  diag$Taxon_simple_name<-record$Taxon_simple_name
  
  # Add slashes to the storage data frame.
  slashes<-rbind(slashes,diag)
  
}

# Backtransform slash y-values.
slashes$y_min<-slashes$y_min^2
slashes$y_max<-slashes$y_max^2

# Create an empty storage data frame for the boxes.
boxes<-data.frame(NULL)

# Loop through each environmental proportional abundance box.
for(i in 1:nrow(props_envi_limits_to_plot)){
  
  # Get the box boundaries.
  record<-props_envi_limits_to_plot[i,]
  
  # Get box sides.
  box_left<-data.frame(x_min=record$x_min,x_max=record$x_min,
                       y_min=record$Quantile_0.025,y_max=record$Quantile_0.975)
  box_right<-data.frame(x_min=record$x_max,x_max=record$x_max,
                        y_min=record$Quantile_0.025,y_max=record$Quantile_0.975)
  box_bottom<-data.frame(x_min=record$x_min,x_max=record$x_max,
                         y_min=record$Quantile_0.025,y_max=record$Quantile_0.025)
  box_top<-data.frame(x_min=record$x_min,x_max=record$x_max,
                      y_min=record$Quantile_0.975,y_max=record$Quantile_0.975)
  
  # Add all box sides to a single data frame.
  box<-rbind(box_left,box_right,box_bottom,box_top)
  
  # Add taxon name to the box.
  box$Taxon_simple_name<-record$Taxon_simple_name
  
  # Add box to the storage box data frame.
  boxes<-rbind(boxes,box)
  
}

# Set proportions plot position dodge value.
props_position_dodge<-0.75

# Plot predicted proportions.
PropsPlot<-ggplot(props_sal,aes(x=Date_numeric,y=Quantile_0.5,color=Life_Stage))+
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
  geom_segment(data=boxes,
               aes(x=x_min,xend=x_max,y=y_min,yend=y_max),
               alpha=0.4,
               color="black",
               linewidth=0.1,
               inherit.aes=F)+
  geom_segment(data=slashes,
               aes(x=x_min,xend=x_max,y=y_min,yend=y_max),
               alpha=0.2,
               color="black",
               linewidth=0.1,
               inherit.aes=F)+
  geom_linerange(aes(ymin=Quantile_0.025,ymax=Quantile_0.975),
                 linewidth=0.5,position=position_dodge(width=props_position_dodge))+
  geom_linerange(aes(ymin=Quantile_0.25,ymax=Quantile_0.75),
                 linewidth=0.75,position=position_dodge(width=props_position_dodge))+
  geom_point(size=1,position=position_dodge(width=props_position_dodge))+
  ggtitle(expression(bold("Proportional Abundances of")~bolditalic("Bd")*bold("-Inhibition Categories on Salamander Skin")),
          subtitle="Point: Posterior median \u2013 Thick line: 50% credible interval \u2013 Thin line: 95% credible interval\nHatched range: Environmental proportional abundance")+
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
        #strip.text=element_text(face=c("italic","plain","plain","plain")))+
  scale_color_manual(values=
                       c("gold","gold3","green2","green4",lighten("dodgerblue1",0.25775,method="absolute",space="HLS"),"dodgerblue3"),
                     name="Stage Class")+
  scale_fill_manual(values=c("white",lighten("gray53",0.25)),name="Site")+
  guides(fill=guide_legend(override.aes=list(color="black",size=0.25)))+
  scale_x_continuous(breaks=unique(props_sal$Date_numeric),
                     labels=levels(props_sal$Date),
                     expand=expansion(mult=c(0,0)))+
  scale_y_sqrt(limits=c(0,NA),expand=expansion(mult=c(0,0.05)))

# Print plot.
print(PropsPlot)

# Save proportions plot.
ggsave(filename="PropsPlot_Bd7.pdf",device=cairo_pdf,plot=PropsPlot,width=14,height=10)

# Done!
print("Done!")
