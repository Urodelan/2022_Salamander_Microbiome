############################
### Check Classification ###
############################

# This script saves plots as pdfs.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load microbe functions.
source("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/HelperFunctions/MicrobeFunctions.R")

# Load sample metadata.
meta<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Microbes/Metadata/SampleMetadata2.csv",stringsAsFactors=F)

# Get names of larval samples.
larval_samples<-meta$SampleID[meta$Type=="Salamander" & meta$LM=="Larvae/Neotene"]

# Get names of metamorphosed samples.
metamorphosed_samples<-meta$SampleID[meta$Type=="Salamander" & meta$LM=="Metamorphosed"]

# Create function for removing tailing numbers from ITS taxa names.
remove_tailing_numbers<-function(x){
  # Split string by semi-colons.
  split<-strsplit(x,split=";")[[1]]
  # Remove tailing number.
  split<-split[-length(split)]
  # Combined string parts separated by semi-colons.
  pasted<-paste0(split,collapse=";")
  # Return new string.
  return(pasted)
}

# Define function for getting taxonomic levels.
getTaxonomicLevel<-function(name,level){
  # Create data frame with translation instructions from level to index.
  taxonomicLevels_translation<-data.frame(
    Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    Index=1:7)
  # Check that the user-specified taxonomic level matches one of the
  # levels in the translation data frame.
  if(!(level %in% taxonomicLevels_translation$Level)) stop(paste0("The specified taxonomic level must be one of the following: ",paste(taxonomicLevels_translation$Level,collapse=", "),"."))
  # Get index of the desired taxonomic level.
  index<-taxonomicLevels_translation$Index[taxonomicLevels_translation$Level==level]
  # Split taxon name by semi-colon.
  x<-strsplit(name,split=";")[[1]]
  # If the taxon is bacterial, remove the trailing number part from the name.
  if(x[1]=="Bacteria") x<-x[-length(x)]
  # If the desired taxonomic level is not included in the taxon name.
  if(length(x) < index){
    # Then assign an NA.
    x<-NA
  } else if(index==7){
    # If the desired taxonomic level is species.
    # Include the genus in the returned taxonomic level as well.
    x<-paste(x[index-1],x[index])
  } else {
    # Otherwise, get the taxonomic level.
    x<-x[index]
  }
  # Return the taxonomic level.
  return(x)
}

# Create storage data frame.
df<-data.frame(NULL)

# Loop through each sample type.
for(j in c("16S","ITS")){
  
  # Read in taxa table.
  reads<-read.taxa.table(file=paste0("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Microbes/Composition/Salamander/Salamander_ASVTable_",j,"_Comp.csv"))
  
  # If the table is for salamander 16S.
  if(j=="16S"){
    # Read in salamander samples which look like negative controls.
    sal_NC<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/ASVTables/Output_16S_2/NC_sal_samples.csv")$NC_sal_samples
    # Remove salamander samples which look like negative controls from the salamander 16S data
    reads<-reads[!(row.names(reads) %in% sal_NC),]
    # Remove fields for salamander 16S taxa which were removed from the dataset with the removal
    # of the samples which look like negative controls.
    reads<-reads[,colSums(reads)!=0]
  }
  
  # If the table is for ITS.
  if(j=="ITS"){
    # Transpose ITS taxa counts.
    transposed<-as.data.frame(t(reads))
    # Create taxon field without tailing numbers.
    transposed$Taxon<-sapply(X=row.names(transposed),FUN=remove_tailing_numbers,USE.NAMES=F)
    # Collapse taxa at highest resolution.
    transposed_reduced<-aggregate(.~Taxon,data=transposed,FUN=sum)
    # Turn taxa names into row names.
    row.names(transposed_reduced)<-transposed_reduced$Taxon
    # Remove taxon as a field.
    transposed_reduced<-transposed_reduced[,-1]
    # Tranpose back to having samples as records.
    reads<-as.data.frame(t(transposed_reduced))
  }
  
  # For larval and metamorphosed samples.
  for(l in c("Larvae/Neotene","Metamorphosed")){
    
    # Get a data frame of read counts for the appropriate samples.
    if(l=="Larvae/Neotene"){
      sub_reads<-reads[row.names(reads) %in% larval_samples,]
    } else {
      sub_reads<-reads[row.names(reads) %in% metamorphosed_samples,]
    }
    
    # Remove taxa which were not observed.
    sub_reads<-sub_reads[,colSums(sub_reads)!=0]
    
    # Get the total read count for each sample.
    total.reads<-rowSums(sub_reads)
    
    # Loop through each taxonomic level below Kingdom.
    for(k in c("Phylum","Class","Order","Family","Genus","Species")){
      # Get the taxonomic level for each taxon.
      taxa<-sapply(X=colnames(sub_reads),FUN=getTaxonomicLevel,level=k,USE.NAMES=F)
      # Subset reads to just those which were classified at the taxonomic level.
      sub<-subset(sub_reads,select=!is.na(taxa))
      # Get the total read count for classified taxa for each sample.
      classified.reads<-rowSums(sub)
      # Get the proportion of classified reads at the taxonomic level for each sample.
      prop.reads.classified<-classified.reads/total.reads
      # Create a data frame with the desired information.
      df_sub<-data.frame(Region=paste(j,"Region"),
                         Life_Stage=l,
                         Level=k,
                         SampleID=row.names(sub_reads),
                         Prop.Classified=prop.reads.classified,
                         row.names=1:nrow(sub_reads))
      # Append the data frame with new information to the storage data frame.
      df<-rbind(df,df_sub)
      # Print progress.
      print(paste0("Completed: ",j,", ",l,", ",k))
    }
    
  }
  
}

# Write out the proportion data frame.
write.csv(df,file="ClassificationProps.csv",row.names=F)

# Read in the proportion data frame.
df<-read.csv(file="ClassificationProps.csv",stringsAsFactors=F)
df$Level<-factor(df$Level,levels=unique(df$Level)) # Order level factor fields.

# Load ggplot.
library(ggplot2)

# Set plot position dodge value.
position_dodge<-0.875

# Plot proportions classified.
Props_plot<-ggplot(df,aes(x=Level,y=Prop.Classified,color=Life_Stage))+
  facet_wrap(facets="Region",
             ncol=1,
             scales="free_y",
             strip.position="top")+
  geom_boxplot(position=position_dodge(width=position_dodge),
               coef=1e6,fill=NA)+
  geom_point(size=1,position=position_dodge(width=position_dodge),alpha=0.5,
             shape=16)+
  ggtitle("Proportion of Salamander Sample Reads Classifed to Taxonomic Level")+
  ylab("Proportion of Reads Classified")+
  xlab(NULL)+
  theme_light()+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.margin=margin(t=-0.1,unit="in"))+
  scale_color_manual(values=c("blue","green4"))+
  labs(color="Life Stage")+
  scale_x_discrete(expand=expansion(mult=c(0.115,0.115)))+
  scale_y_continuous(expand=expansion(mult=c(0.025,0.025)),
                     limits=c(0,1))

# Print plot.
print(Props_plot)

# Save plot.
ggsave(filename="Classified_plot.pdf",plot=Props_plot,width=8,height=7,units="in")

# Add a site field to the data.
df$Site<-meta$Site[match(df$SampleID,meta$SampleID)]

# Set plot position dodge value.
position_dodge2<-0.875

# Plot proportions classified.
Props_plot2<-ggplot(df,aes(x=Level,y=Prop.Classified,color=Life_Stage,shape=Site,linetype=Site))+
  facet_wrap(facets="Region",
             ncol=1,
             scales="free_y",
             strip.position="top")+
  geom_boxplot(position=position_dodge(width=position_dodge2),
               coef=1e6,fill=NA)+
  geom_point(size=1,position=position_dodge(width=position_dodge2),alpha=0.5)+
  ggtitle("Proportion of Reads in Each Salamander Sample Classifed Down to At Least the Given Taxonomic Level")+
  ylab("Proportion of Reads Classified")+
  xlab(NULL)+
  theme_light()+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="bottom",
        legend.margin=margin(t=-0.05,l=0.375,r=0.375,unit="in"))+
  scale_color_manual(values=c("blue","green4"))+
  labs(linetype="Site:",shape="Site:",color="Life Stage:")+
  scale_x_discrete(expand=expansion(mult=c(0.115,0.115)))+
  scale_y_continuous(expand=expansion(mult=c(0.025,0.025)),
                     limits=c(0,1))+
  scale_shape_manual(values=c(16,8))

# Print plot.
print(Props_plot2)

# Save plot.
ggsave(filename="Classified_plot2.pdf",plot=Props_plot2,width=10,height=8,units="in")

# Done!
print("Done!")
