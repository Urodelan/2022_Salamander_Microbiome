# Use PCAs resampling salamanders to 4 random subsets samples (and using proportional abundances)
# to identify salamander samples which are similar to negative controls.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###############################
### Check ASV Table for 16S ###
###############################

# Set seed.
set.seed(1234)

# Load 16S ASV table.
tb16S<-read.csv("../ASV_Table_16S.csv",row.names=1)

# Define function for simplexing a table.
simplex_table<-function(table){
  simplexed_table<-as.data.frame(t(apply(X=table,MARGIN=1,FUN=function(x) x/sum(x))))
  return(simplexed_table)
}

# Add a sample field to the 16S ASV table.
tb16S$SampleID<-paste0(sapply(strsplit(row.names(tb16S),split="_"),"[[",1),"_",sapply(strsplit(row.names(tb16S),split="_"),"[[",2))

# Collapse 16S ASV table by sample.
tb16S_agg<-aggregate(.~SampleID,data=tb16S,FUN=sum)

# Load sample metadata.
meta<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Metadata/SampleMetadata.csv")

# Remove the samples P5_13, P3_11, G6_23, and CB_10 from the data.
tb16S_agg<-tb16S_agg[!(tb16S_agg$SampleID %in% c("P5_13","P3_11","G6_23","CB_10")),]

# Remove the samples P5_13, P3_11, G6_23, and CB_10 from the metadata.
meta<-meta[!(meta$SampleID %in% c("P5_13","P3_11","G6_23","CB_10")),]

# Get 16S ASV table with just biological sequences.
tb16S_bio<-tb16S_agg[,!(grepl("^Coligo",colnames(tb16S_agg)) | colnames(tb16S_agg)=="Synthgene16S")]
row.names(tb16S_bio)<-tb16S_bio$SampleID
tb16S_bio<-tb16S_bio[,-which(colnames(tb16S_bio)=="SampleID")]

# Turn 16S ASV table into proportional abundances (excluding the synthgene).
tb16S_bio_props<-simplex_table(tb16S_bio)

# Reorder proportional abundances (excluding synthgene) to be the same as the metadata.
tb16S_bio_props<-tb16S_bio_props[match(meta$SampleID,row.names(tb16S_bio_props)),]

# Check that tables and metadata are lined up for proportional abundance data.
identical(row.names(tb16S_bio_props),as.character(meta$SampleID)) # Yes.

# Change up water names to indicate whether the sample is a control.
meta$Type<-as.character(meta$Type)
meta$Type<-ifelse(meta$Type=="Water" & meta$Control=="Yes","Water_Control",meta$Type)
meta$Type<-as.factor(meta$Type)

# Replace sample type underscores with spaces.
meta$Type<-gsub(pattern="_",replacement=" ",x=meta$Type)

# Reorder sample type levels.
meta$Type<-factor(meta$Type,levels=c("Salamander","Wet Swab","Dry Swab","Blank","Water","Water Control","Substrate","Mock Community"))

# Get summary of sample type counts.
print(table(meta$Type))

# Design four random subsets for the salamander samples.
sal_rand_subs<-rep(1:4,length.out=sum(meta$Type=="Salamander"))
sal_rand_subs<-sal_rand_subs[sample(x=1:length(sal_rand_subs),size=length(sal_rand_subs),replace=F)]

# Get row indices of salamander samples.
sal_indices<-which(meta$Type=="Salamander")

# Subset metadata to just salamander samples.
meta_sal<-meta[sal_indices,]

# Subset proportional abundances to just salamander samples.
sal_prop<-tb16S_bio_props[sal_indices,]

# Subset metadata to just non-salamander samples.
meta_nonsal<-meta[-sal_indices,]

# Subset proportional abundances to just non-salamander samples.
nonsal_prop<-tb16S_bio_props[-sal_indices,]

# Create an empty list for storing salamander subsample plots.
plot_list<-vector(mode="list",length=4)
plot_list2<-vector(mode="list",length=4)

# Create and store 4 subsample plots.

# Create data frame specifying PC1 and PC2 coordinate boundaries for
# each subsetted PCA for identifying salamander samples which are
# similar to negative controls.
NC_coord_bounds<-data.frame(Sub=1:4,
                            PC1_bound_val=c(-0.05,-0.05,0,0),
                            PC1_bound_type=c("max","max","max","min"),
                            PC2_bound_val=c(-0.2,-0.2,0.2,0.25),
                            PC2_bound_type=c("max","max","min","min"))

# Initiate empty storage vector for holding names of salamander samples
# which look like negative control samples on the PCA subsets.
NC_sal_samples<-c()

## Loop through 4 iterations.
library(ggbiplot)
for(i in 1:4){
  
  # Get indices for salamander sample subset.
  sal_sub_indices<-which(sal_rand_subs==i)

  # Subset salamander sample metadata to the sample subset.
  meta_sal_sub<-meta_sal[sal_sub_indices,]
  
  # Subset salamander sample proportional abundances to the sample subset.
  sal_prop_sub<-sal_prop[sal_sub_indices,]
  
  # Combine non-salamander and salamander sample subset.
  full_sub_meta<-rbind(meta_sal_sub,meta_nonsal)
  
  # Combine non-salamander and salamander sample subset proportional abundances.
  full_sub_props<-rbind(sal_prop_sub,nonsal_prop)
  
  # Visualize salamander subsets sample type with PCA with labels.
  PCA<-prcomp(full_sub_props)
  g<-ggbiplot(PCA,obs.scale=1,var.scale=1,groups=full_sub_meta$Type,
              ellipse=TRUE,ellipse.prob=0.95,alpha=0.5,labels=row.names(full_sub_props))
  g<-g+theme_light()+
    # scale_color_discrete(name='Sample Type')+
    ggtitle("PCA on Sample Type",
            subtitle=paste0("16S Region - Subsample #",i))+
    theme(legend.direction="horizontal",legend.position="bottom",
          plot.title=element_text(hjust=0.5,face="bold"),
          plot.subtitle=element_text(hjust=0.5))+
    guides(color=guide_legend(nrow=2,byrow=T))+
    scale_color_manual(name='Sample Type',
                       values=c("green4","gold","gold3","purple","blue","dodgerblue1","brown","pink2"))
  g$layers<-c(g$layers[[2]],g$layers[[3]])
  
  # Visualize salamander subsets sample type with PCA without labels.
  g2<-ggbiplot(PCA,obs.scale=1,var.scale=1,groups=full_sub_meta$Type,
              ellipse=TRUE,ellipse.prob=0.95,alpha=0.5)
  g2<-g2+theme_light()+
    # scale_color_discrete(name='Sample Type')+
    ggtitle("PCA on Sample Type",
            subtitle=paste0("16S Region - Subsample #",i))+
    theme(legend.direction="horizontal",legend.position="bottom",
          plot.title=element_text(hjust=0.5,face="bold"),
          plot.subtitle=element_text(hjust=0.5))+
    guides(color=guide_legend(nrow=2,byrow=T))+
    scale_color_manual(name='Sample Type',
                       values=c("green4","gold","gold3","purple","blue","dodgerblue1","brown","pink2"))
  g2$layers<-c(g2$layers[[2]],g2$layers[[3]])
  
  # Get PC1 and PC2 from the PCA.
  coords<-as.data.frame(PCA$x[,1:2])
  
  # Add sample type to the PCA coordinates.
  coords$Type<-meta$Type[match(row.names(coords),meta$SampleID)]
  
  # Subset PCA coordinates to just salamander samples.
  coords<-coords[coords$Type=="Salamander",]
  
  # Subset the negative control PCA bounding coordinates to the current PCA subset.
  NC_coord_bounds_sub<-NC_coord_bounds[NC_coord_bounds$Sub==i,]
  
  # Get just PCA points which meet the PC1 bounding criterion.
  if(NC_coord_bounds_sub$PC1_bound_type=="min"){
    coords<-coords[coords$PC1 > NC_coord_bounds_sub$PC1_bound_val,]
  } else {
    coords<-coords[coords$PC1 < NC_coord_bounds_sub$PC1_bound_val,]
  }
  
  # Get just PCA points which meet the PC2 bounding criterion.
  if(NC_coord_bounds_sub$PC2_bound_type=="min"){
    coords<-coords[coords$PC2 > NC_coord_bounds_sub$PC2_bound_val,]
  } else {
    coords<-coords[coords$PC2 < NC_coord_bounds_sub$PC2_bound_val,]
  }
  
  # Store names of salamander samples which look like negative controls on the PCA subset.
  NC_sal_samples<-c(NC_sal_samples,row.names(coords))
  
  # Store labeled PCA plot in the storage list.
  plot_list[[i]]<-g
  
  # Store unlabeled PCA plot in the storage list.
  plot_list2[[i]]<-g2
  
}

# Arrange plots in a panel.
library(gridExtra)
plot_panel<-arrangeGrob(grobs=plot_list,nrow=2,ncol=2)
plot_panel2<-arrangeGrob(grobs=plot_list2,nrow=2,ncol=2)

# Print plot panel.
grid.arrange(plot_panel)
grid.arrange(plot_panel2)

# Save PCA plot panel.
ggsave(filename="PCA_SampleType_16S_Subsamples_Labeled.tiff",plot=plot_panel,width=14,height=10,units="in",dpi=300)
ggsave(filename="PCA_SampleType_16S_Subsamples2.tiff",plot=plot_panel2,width=14,height=10,units="in",dpi=300)

# Write out salamander samples which look like negative controls on the subsetted PCA plots.
write.csv(as.data.frame(NC_sal_samples),file="NC_sal_samples.csv",row.names=F)

# Done!
print("Done!")
