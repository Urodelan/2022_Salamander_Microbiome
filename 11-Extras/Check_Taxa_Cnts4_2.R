#########################
### Check Taxa Counts ###
#########################

# This script produces an updated plot from Check_Taxa_Cnts3.R.

# This script saves plots as pdfs.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load microbe functions.
source("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/HelperFunctions/MicrobeFunctions.R")

# Create storage data frame.
df<-data.frame(NULL)

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

# Loop through each sample type.
for(k in c("Salamander","Substrate","Water")){
  # Loop through each barcode region.
  for(j in c("16S","ITS")){
    # Read in taxa table.
    reads<-read.taxa.table(file=paste0("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/Microbes/Composition/",k,"/",k,"_ASVTable_",j,"_Comp.csv"))
    # If the table is for salamander 16S.
    if(k=="Salamander" & j=="16S"){
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
    # Rename columns.
    colnames(reads)<-paste0("T",1:ncol(reads))
    # Loop through each record and sort values in descending order.
    for(i in 1:nrow(reads)){
      reads[i,]<-sort(reads[i,],decreasing=T)
    }
    # Copy reads.
    cum_reads<-reads
    # Loop through each column.
    for(i in 2:ncol(cum_reads)){
      # Get cumulative read counts.
      cum_reads[,i]<-rowSums(reads[,1:i])
    }
    # Loop through each record.
    for(i in 1:nrow(cum_reads)){
      # Loop through each column.
      for(m in 1:ncol(cum_reads)){
        # Calculate cumulative proportion of reads.
        cum_reads[i,m]<-cum_reads[i,m]/cum_reads[i,ncol(cum_reads)]
      }
    }
    # Calculate the median and IQR of cumulative proportions.
    cum_props_metrics<-as.data.frame(t(apply(X=cum_reads,MARGIN=2,FUN=quantile,probs=c(0.25,0.5,0.75))))
    ## Rename columns.
    colnames(cum_props_metrics)<-c("Q1","Median","Q3")
    ## Add field for cumulative number of taxa.
    cum_props_metrics$Cum_Num_Taxa<-1:nrow(cum_props_metrics)
    ## Add in a record for zero cumulative taxa.
    cum_props_metrics<-rbind(data.frame(Q1=0,Median=0,Q3=0,Cum_Num_Taxa=0),
                             cum_props_metrics)
    # Add field denoting sample type.
    cum_props_metrics$Type<-k
    # Add field denoting barcode region.
    cum_props_metrics$Region<-j
    # Re-order fields.
    cum_props_metrics<-cum_props_metrics[,c("Region","Type","Cum_Num_Taxa","Q1","Median","Q3")]
    # Rename rows.
    row.names(cum_props_metrics)<-1:nrow(cum_props_metrics)
    # Add data to the storage data frame.
    df<-rbind(df,cum_props_metrics)
    # Print progress.
    print(paste(k,j,"completed."))
  }
}

# Write out cumulative curve data frame.
write.csv(df,file="Cum_curve3.csv",row.names=F)

# Read in cumulative curve data frame.
df<-read.csv("Cum_curve3.csv")

# Remove records where the cumulative number of taxa is zero.
df<-df[df$Cum_Num_Taxa!=0,]

# Create factor levels for the cumulative curve data frame.
df$Type<-factor(df$Type,levels=c("Salamander","Water","Substrate"))

# Add region to the end of 16S and ITS.
df$Region<-paste(df$Region,"Region")

# Load ggplot.
library(ggplot2)

# Create major and minor breaks between 1 and 100,000.
plot_major_breaks<-10^seq(0,5,by=1)
plot_minor_breaks<-c()
for(i in 1:(length(plot_major_breaks)-1)){
  plot_minor_breaks<-c(plot_minor_breaks,seq(plot_major_breaks[i],plot_major_breaks[i+1],by=plot_major_breaks[i]))
}
plot_minor_breaks<-plot_minor_breaks[!duplicated(plot_minor_breaks)]

# Create cumulative read plot #3.
cum_read_plot_4<-ggplot(data=df,aes(x=Cum_Num_Taxa,color=Type))+
  facet_wrap(facets="Region",
             ncol=2,
             scales="free_x",
             strip.position="top")+
  geom_line(aes(y=Q1),alpha=0.75,linetype="dotted")+
  geom_line(aes(y=Q3),alpha=0.75,linetype="dotted")+
  geom_line(aes(y=Median),alpha=0.75)+
  scale_color_manual(values=c("green4","blue","brown"))+
  theme_light()+
  ggtitle("Cumulative Number of Taxa and Proportion of Reads",
          subtitle="Solid line: Median          Dotted lines: Interquartile range")+
  xlab("Cumulative Number of Sorted Taxa")+
  ylab("Cumulative Proportion of Reads")+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.margin=margin(t=-0.1,unit="in"),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1))+
  scale_y_continuous(expand=expansion(mult=c(0,0.025)),
                     limits=c(0,1))+
  scale_x_log10(expand=expansion(mult=c(0,0)), # Note the log10 x-axis scale.
                limits=c(1,NA),
                breaks=plot_major_breaks,
                minor_breaks=plot_minor_breaks)+
  # scale_x_continuous(expand=expansion(mult=c(0.01,0)),
  #                    limits=c(0,NA))+
  guides(color=guide_legend(order=1),linetype=guide_legend(order=2))

# Display plot.
plot(cum_read_plot_4)

# Save plot.
ggsave(filename="Cum_reads_plot_5.pdf",plot=cum_read_plot_4,
       width=7,height=5,units="in")

# Done!
print("Done!")
