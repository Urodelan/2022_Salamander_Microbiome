################################
### Plot ITS Summary Indices ###
################################

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

# Read in Hill's diversity predictions.
hills<-read.csv(file="../FromCluster/Predictions/hills2_summary_ITS.csv")

# Format date field.
hills$Date<-as.Date(hills$Date)

# Re-format the date field.
hills$Date<-format(hills$Date,"%b %e")

# Remove extra space from the date field.
hills$Date<-gsub("  "," ",hills$Date)

# Format date field as a factor.
hills$Date<-factor(hills$Date,levels=unique(hills$Date))

# Create a numeric date field.
hills$Date_numeric<-as.numeric(hills$Date)

# Create life stage field.
hills$Life_Stage<-paste(hills$Age,ifelse(hills$LM=="Larvae/Neotene",ifelse(hills$Age=="Age-2+","Neotene","Larvae"),"Metamorphosed"))

# Format life stage field as a factor.
hills$Life_Stage<-factor(hills$Life_Stage,levels=c("Age-0 Larvae","Age-0 Metamorphosed","Age-1 Larvae","Age-1 Metamorphosed","Age-2+ Neotene","Age-2+ Metamorphosed"))

# Create data frame storing information on shading by site.
shading<-data.frame(min=seq(from=0.5,to=max(hills$Date_numeric)-0.5,by=1),
                    max=seq(from=1.5,to=max(hills$Date_numeric)+0.5,by=1),
                    col=c("Gibson Lakes","Ponds Lake"))

# Load ggplot2.
library(ggplot2)

# Load package containing lighten function.
library(colorspace)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Set Hill's plot position dodge value.
hills_position_dodge<-0.75

# Plot Hill's diversity index.
hillsPlot<-ggplot(hills,aes(x=Date_numeric,y=Quantile_0.5,color=Life_Stage))+
  geom_rect(data=shading,
            aes(xmin=min,xmax=max,ymin=-Inf,ymax=Inf,fill=factor(col)),
            alpha=0.4,
            inherit.aes=F)+
  geom_linerange(aes(ymin=Quantile_0.025,ymax=Quantile_0.975),
                 size=0.5,position=position_dodge(width=hills_position_dodge))+
  geom_linerange(aes(ymin=Quantile_0.25,ymax=Quantile_0.75),
                 size=1,position=position_dodge(width=hills_position_dodge))+
  geom_point(size=1.5,position=position_dodge(width=hills_position_dodge))+
  ggtitle("Fungal Hill's Diversity Index (\u03b1 = 2) on Salamander Skin",
          subtitle="Point: Posterior median\nThick line: 50% credible interval \u2013 Thin line: 95% credible interval")+
  ylab("Hill's Diversity Index")+
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
  scale_x_continuous(breaks=unique(hills$Date_numeric),
                     labels=levels(hills$Date),
                     expand=expansion(mult=c(0,0)))

# Display Hill's diversity plot.
print(hillsPlot)

# Save Hill's diversity plot.
ggsave(filename="HillsPlot_ITS.pdf",device=cairo_pdf,plot=hillsPlot,width=7.5,height=5)

# Load Bray-Curtis dissimilarity.
load("../FromCluster/Predictions/bc_summary_ITS.RData")

# Load additional plotting packages.
library(lattice)
library(cowplot)
library(gridGraphics)
library(plot3D)

# Create empty storage list for Bray-Curtis summary plots for each life stage.
plot_rows<-vector(mode="list",length(bc_summary))

# Loop through each stage class.
for(i in 1:length(bc_summary)){
  # Get Bray-Curtis summaries for the stage class.
  bc_sub<-bc_summary[i]
  # Get the stage class name.
  bc_stage_class<-names(bc_sub)
  # Get the minimum BC value for the stage class, rounded down to the nearest hundredth.
  scale_min<-floor(min(unlist(bc_sub),na.rm=T)*100)/100
  # Get the maximum BC value for the stage class, rounded up to the nearest hundredth.
  scale_max<-ceiling(max(unlist(bc_sub),na.rm=T)*100)/100
  # Get the dates.
  bc_dates<-colnames(bc_sub[[1]][[1]])
  # Get the sites associated with the dates.
  bc_dates_sites<-sapply(X=strsplit(x=bc_dates,split="_"),FUN="[[",1)
  # Assign date colors based on the site.
  date.colors<-ifelse(bc_dates_sites=="Gibson","blue","green4")
  # Get the dates without site information.
  bc_dates<-sapply(X=strsplit(x=bc_dates,split="_"),FUN="[[",2)
  # Format the dates.
  bc_dates<-as.Date(bc_dates)
  # Re-format the dates.
  bc_dates<-format(bc_dates,"%b %e")
  # Remove extra space from the dates.
  bc_dates<-gsub("  "," ",bc_dates)
  # Create a storage list for BC quantile plots for the stage class.
  bc_sub_plots<-vector(mode="list",length(bc_sub[[1]])+2)
  # Loop through each quantile.
  for(j in 1:length(bc_sub[[1]])){
    # Get the quantile.
    bc_quant<-as.matrix(bc_sub[[1]][[j]])
    # Get the quantile name.
    bc_quant_name<-as.numeric(strsplit(x=names(bc_sub[[1]])[j],split="_")[[1]][2])*100
    # Plot the quantile.
    bc_sub_plots[[j+1]]<-levelplot(bc_quant,xlab=NULL,ylab=NULL,
              #main=list(label=paste0("     ",bc_quant_name," percentile"),fontface="plain",cex=1),
              col.regions=heat.colors(1000,rev=T),
              ylim=c(ncol(bc_quant)+0.5,0.5),xlim=c(0.5,ncol(bc_quant)+0.5),
              scale=list(
                x=list(
                  alternating=2,
                  col=date.colors,
                  at=c(1:ncol(bc_quant)),
                  labels=bc_dates,
                  cex=0.5,rot=90),
                y=list(
                  col=date.colors,
                  at=c(1:ncol(bc_quant)),
                  labels=bc_dates,
                  cex=0.5)
              ),
              at=seq(scale_min,scale_max,by=0.001),
              colorkey=F,
              par.settings=list(panel.background=list(col=alpha("lightgrey",0.35))))
  }
  
  # Create a color key legend.
  ## Create a large right margin (last digit).
  par(mar=c(5,0,5,12))
  ## Create a color key which is wide relative to its height.
  colkey(col=heat.colors(1000,rev=T),clim=c(scale_min,scale_max),length=1,width=5)
  ## Record the plot in the storage list.
  bc_sub_plots[[length(bc_sub_plots)]]<-recordPlot()
  ## Reset margins back to default settings.
  par(mar=c(5.1,4.1,4.1,2.1))
  
  # Make a graphical label object.
  bc_sub_plots[[1]]<-ggdraw()+
    draw_label(bc_stage_class,hjust=0.5,x=0.5,y=0.5,angle=90)+
    theme(text=element_text(family="Arial"))

  # Put all Bray-Curtis dissimilarity plots together.
  # Use smaller relative heights and widths for the color key legend.
  plot_rows[[i]]<-plot_grid(plotlist=bc_sub_plots,ncol=length(bc_sub_plots),
                            rel_widths=c(0.2,rep(1,5),0.5),
                            rel_heights=c(0.5,rep(1,5),0.5))
  
}

# Combine BC plots from all stage classes.
BC_all_stages<-plot_grid(plotlist=plot_rows,nrow=length(plot_rows))

# Set number of spaces between percentile text strings.
BC_num_spaces<-24

# Set number of spaces to append to the end of the overall text string.
# The higher this value, the further to the left the overall text string is pushed.
BC_num_spaces_at_end<-7

# Make a graphical supertitle object.
BC_title<-ggdraw()+
  draw_label("Fungal Bray-Curtis Dissimilarity Through Time on Salamander Skin",fontface="bold",
             hjust=0.5,x=0.5,y=0.5)+
  draw_label("Blue dates: Gibson Lakes \u2013 Green dates: Ponds Lake",x=0.5,y=0.25)+
  draw_label(paste0("2.5 percentile",paste0(rep(" ",BC_num_spaces),collapse=""),
                    "25 percentile",paste0(rep(" ",BC_num_spaces),collapse=""),
                    "50 percentile",paste0(rep(" ",BC_num_spaces),collapse=""),
                    "75 percentile",paste0(rep(" ",BC_num_spaces),collapse=""),
                    "97.5 percentile",paste0(rep(" ",BC_num_spaces_at_end),collapse="")),
             x=0.5,y=0)+
  theme(text=element_text(family="Arial"))

# Add supertitle to the plot panel.
plot_BC<-plot_grid(BC_title,BC_all_stages,ncol=1,rel_heights=c(0.1,1))

# Save Bray-Curtis dissimilarity index plot.
ggsave("BrayCurtisPlot_ITS.pdf",device=cairo_pdf,plot_BC,width=14,height=15)

# Done!
print("Done!")
