# Re-doing the PCAs with proportional abundances instead of raw counts.

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###############################
### Check ASV Table for 16S ###
###############################

# Load 16S ASV table.
tb16S<-read.csv("ASV_Table_16S.csv",row.names=1)

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

# Are there any samples which are not included in the metadata?
sum(!(tb16S_agg$SampleID %in% meta$SampleID)) # Nope.

# Are there any field samples that did not amplify?
sum(!(meta$SampleID %in% tb16S_agg$SampleID)) # Nope.

# Get 16S ASV table with only non-biological sequences.
tb16S_nonbio<-tb16S_agg[,colnames(tb16S_agg)=="SampleID" | grepl("^Coligo",colnames(tb16S_agg)) | grepl("^Synthgene",colnames(tb16S_agg))]

# Read in sample positions.
samplePositions16S<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Lab/LibraryPrep/BarcodePairCSVs/GompertPrep1_03-03-2020_16S.csv")

# Read in 16S coligos with corrected positions.
coligos_16S_corr<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Lab/LibraryPrep/Coligos16S_ProperPositions.csv")

# Read in 16S coligos with uncorrected positions.
coligos_16S_uncorr<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Lab/LibraryPrep/Coligos16S.csv")

# Check that coligos are the same between files.
identical(coligos_16S_uncorr$SequenceWithAmpPrimers,coligos_16S_corr$SequenceWithAmpPrimers) # True.

# Create link between incorrect old coligo names and new correct coligo names.
coligos_16S<-data.frame(Old_Name=paste0("Coligo16S_",coligos_16S_uncorr$Position),New_Name=paste0("Coligo16S_",coligos_16S_corr$Position))

# Replace old incorrect coligo names with new correct coligo names.
for(i in 2:(ncol(tb16S_nonbio)-1)){
  colnames(tb16S_nonbio)[i]<-as.character(coligos_16S$New_Name[coligos_16S$Old_Name==colnames(tb16S_nonbio)[i]])
}

# Check if all samples have a single well position.
all(aggregate(Well.Position~sample.name,data=samplePositions16S,FUN=function(x) length(unique(x)))$Well.Position==1) # False.

# What are the number of well positions samples have?
unique(aggregate(Well.Position~sample.name,data=samplePositions16S,FUN=function(x) length(unique(x)))$Well.Position) # 1 and 2.

# Which samples have more than one well position?
aggregate(Well.Position~sample.name,data=samplePositions16S,FUN=function(x) length(unique(x)))$sample.name[aggregate(Well.Position~sample.name,data=samplePositions16S,FUN=function(x) length(unique(x)))$Well.Position==2] # Samples from plate KG5.

# View the well positions of the samples from plate KG5.
samplePositions16S[samplePositions16S$sample.name %in% c("G8_1S","G8_2S","G8_3S","G8_3W","G8_4S","G8_4W","MC_1","MC_2"),]

# Check for instances where coligos from improper cell appear in the samples.
cell_contam_sample_name<-c()
cell_contam_position<-c()
cell_contam_amount<-c()
cell_contam_proper<-c()
for(i in 1:nrow(tb16S_nonbio)){
  sample_name<-tb16S_nonbio$SampleID[i]
  sample_positions<-as.character(unique(samplePositions16S$Well.Position[samplePositions16S$sample.name==sample_name]))
  for(j in 2:(ncol(tb16S_nonbio)-1)){
    if(tb16S_nonbio[i,j]>0){
      position<-strsplit(colnames(tb16S_nonbio)[j],split="_")[[1]][2]
      if(!(position %in% sample_positions)){
        cell_contam_sample_name<-c(cell_contam_sample_name,sample_name)
        cell_contam_position<-c(cell_contam_position,position)
        cell_contam_amount<-c(cell_contam_amount,tb16S_nonbio[i,j])
        cell_contam_proper<-c(cell_contam_proper,sum(tb16S_nonbio[i,colnames(tb16S_nonbio) %in% paste0("Coligo16S_",sample_positions)]))
      }
    }
  }
}
cell_contamination<-data.frame(SampleID=cell_contam_sample_name,ImproperCell=cell_contam_position,Quantity=cell_contam_amount,AmtOfProper=cell_contam_proper)
cell_contamination$Ratio<-cell_contamination$Quantity/cell_contamination$AmtOfProper

# Subset to samples which are of special concern.
checkThese<-subset(cell_contamination,Ratio>0.1) # Be especially wary of samples: P5_13, P3_11, G6_23, and CB_10.

# Check that all samples have biological amplification.
tb16S_for_checking_bio_amp<-tb16S_agg[,!(grepl("^Coligo",colnames(tb16S_agg)) | grepl("^Synthgene",colnames(tb16S_agg)))]
row.names(tb16S_for_checking_bio_amp)<-tb16S_for_checking_bio_amp$SampleID
tb16S_for_checking_bio_amp<-tb16S_for_checking_bio_amp[,-which(colnames(tb16S_for_checking_bio_amp)=="SampleID")]
## All samples have biological amplification.
row.names(tb16S_for_checking_bio_amp)[rowSums(tb16S_for_checking_bio_amp)==0]

# Remove the samples P5_13, P3_11, G6_23, and CB_10 from the data.
tb16S_agg<-tb16S_agg[!(tb16S_agg$SampleID %in% c("P5_13","P3_11","G6_23","CB_10")),]

# Remove the samples P5_13, P3_11, G6_23, and CB_10 from the metadata.
meta<-meta[!(meta$SampleID %in% c("P5_13","P3_11","G6_23","CB_10")),]

# Get 16S ASV table without any coligos sequences.
tb16S_bio<-tb16S_agg[,!grepl("^Coligo",colnames(tb16S_agg))]
row.names(tb16S_bio)<-tb16S_bio$SampleID
tb16S_bio<-tb16S_bio[,-which(colnames(tb16S_bio)=="SampleID")]

# Are there any samples without any synthgene amplified?
sum(tb16S_bio$Synthgene16S==0) # Yes, 13 samples.

# What are the samples without any synthgene amplified?
row.names(tb16S_bio)[tb16S_bio$Synthgene16S==0] # All are substrate and water samples.

# Let's try adding one to these synthgene counts which were not detected.
tb16S_bio$Synthgene16S[tb16S_bio$Synthgene16S==0]<-1

# Turn 16S ASV table into proportional abundances (including the synthgene).
tb16S_bio_props<-as.data.frame(t(apply(tb16S_bio,MARGIN=1,FUN=function(x) x/sum(x))))

# Scale biological sequence abundances by the synthgene abundance.
tb16S_bio_scl<-as.data.frame(t(apply(tb16S_bio_props,MARGIN=1,FUN=function(x) x/x["Synthgene16S"])))

# Remove the synthgene from the dataset.
tb16S_bio_scl<-tb16S_bio_scl[,-which(colnames(tb16S_bio_scl)=="Synthgene16S")]

# Reorder scaled abundances to be the same as the metadata.
tb16S_bio_scl<-tb16S_bio_scl[match(meta$SampleID,row.names(tb16S_bio_scl)),]

# Turn 16S ASV table into proportional abundances (excluding the synthgene).
tb16S_bio_props<-simplex_table(tb16S_bio[,-which(colnames(tb16S_bio)=="Synthgene16S")])

# Reorder proportional abundances (excluding synthgene) to be the same as the metadata.
tb16S_bio_props<-tb16S_bio_props[match(meta$SampleID,row.names(tb16S_bio_props)),]

# Get table of raw sequence counts.
tb16S_bio<-tb16S_bio[,-which(colnames(tb16S_bio)=="Synthgene16S")]

# Reorder raw sequence counts to be the same as the metadata.
tb16S_bio<-tb16S_bio[match(meta$SampleID,row.names(tb16S_bio)),]

# Check that tables and metadata are lined up.
## Raw count data.
identical(row.names(tb16S_bio),as.character(meta$SampleID)) # Yes.
## Proportional abundances data.
identical(row.names(tb16S_bio_props),as.character(meta$SampleID)) # Yes.
## Absolute abundances data.
identical(row.names(tb16S_bio_scl),as.character(meta$SampleID)) # Yes.

# Check that the number of columns are the same between all tables.
ncol(tb16S_bio); ncol(tb16S_bio_props); ncol(tb16S_bio_scl)

# Remove the samples P5_13, P3_11, G6_23, and CB_10 from the un-aggregated data.
tb16S<-tb16S[!(tb16S$SampleID %in% c("P5_13","P3_11","G6_23","CB_10")),]
tb16S<-tb16S[,-which(colnames(tb16S)=="SampleID")]

# Remove non-biological sequences from the un-aggregated data.
tb16S<-tb16S[,!(grepl("^Coligo",colnames(tb16S)) | grepl("^Synthgene",colnames(tb16S)))]

# Remove records without any biological sequences from un-aggregated data.
tb16S<-tb16S[rowSums(tb16S)!=0,]

# Get PCR number from the 16S ASV table.
tb16S_PCR<-as.factor(substr(row.names(tb16S),start=nchar(row.names(tb16S)),stop=nchar(row.names(tb16S))))

# Visualize PCR with PCA.
library(ggbiplot)
PCA<-prcomp(simplex_table(tb16S))
g<-ggbiplot(PCA,obs.scale=1,var.scale=1,groups=tb16S_PCR,
            ellipse=TRUE,ellipse.prob=0.95,alpha=0.5)
g<-g+theme_light()+
  scale_color_discrete(name='PCR Number')+
  ggtitle("PCA on PCR",
          subtitle="16S Region")+
  theme(legend.direction="horizontal",legend.position="bottom",
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5))+
  guides(color=guide_legend(nrow=1))
g$layers<-c(g$layers[[2]],g$layers[[3]])
print(g)

# Save PCA plot.
ggsave(filename="Output_16S_2/PCA_PCR_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# Review of tables...
## Raw counts: tb16S_bio
## Proportional abundances: tb16S_bio_props
## Absolute abundances: tb16S_bio_scl

# Change up water names to indicate whether the sample is a control.
meta$Type<-as.character(meta$Type)
meta$Type<-ifelse(meta$Type=="Water" & meta$Control=="Yes","Water_Control",meta$Type)
meta$Type<-as.factor(meta$Type)

# Replace sample type underscores with spaces.
meta$Type<-gsub(pattern="_",replacement=" ",x=meta$Type)

# Reorder sample type levels.
meta$Type<-factor(meta$Type,levels=c("Salamander","Wet Swab","Dry Swab","Blank","Water","Water Control","Substrate","Mock Community"))

# Visualize sample type with PCA.
library(ggbiplot)
PCA<-prcomp(tb16S_bio_props)
g<-ggbiplot(PCA,obs.scale=1,var.scale=1,groups=meta$Type,
            ellipse=TRUE,ellipse.prob=0.95,alpha=0.5)
g<-g+theme_light()+
  # scale_color_discrete(name='Sample Type')+
  ggtitle("PCA on Sample Type",
          subtitle="16S Region")+
  theme(legend.direction="horizontal",legend.position="bottom",
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5))+
  guides(color=guide_legend(nrow=2,byrow=T))+
  scale_color_manual(name='Sample Type',
                     values=c("green4","gold","gold3","purple","blue","dodgerblue1","brown","pink2"))
g$layers<-c(g$layers[[2]],g$layers[[3]])
print(g)

# Save PCA plot.
ggsave(filename="Output_16S_2/PCA_SampleType_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# Get raw read counts for sample types.
read_quantity_raw<-rowSums(tb16S_bio)
read_quantity_raw<-data.frame(SampleID=names(read_quantity_raw),Quantity=read_quantity_raw)

# Add types to salamander read quantity.
read_quantity_raw$Type<-meta$Type

# Quick plot of raw abundances.
library(ggplot2)
(g<-ggplot(data=read_quantity_raw,aes(x=Type,y=Quantity))+
    geom_violin(bw=0.375,trim=F)+
    geom_point(pch=21,alpha=0.5)+
    ylab("Read Count")+
    ggtitle("Read Count by Sample Type",
            subtitle="16S Region")+
    theme_light()+
    xlab("Sample Type")+
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1),
          plot.subtitle=element_text(hjust=0.5))+
    scale_y_log10())

# Save raw abundance plot.
ggsave(filename="Output_16S_2/ReadCount_SampleType_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# From full lab prep protocol:
## 3 ul of 0.03 pg/ul = 0.09 pg synthgene added to each sample.

# What is the molecular weight of the synthgenes?

## With attached amplification primers.
library(Biostrings)
(synthgene16S_with_amp<-paste0("GTGCCAGCAGCCGCGGTAA","GCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT",as.character(reverseComplement(DNAString("GGACTACTAGGGTATCTAAT")))))
(synthgeneITS_with_amp<-paste0("CTTGGTCATTTAGAGGAAGTAA","TGCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT",as.character(reverseComplement(DNAString("GCTGCGTTCTTCATCGATGC")))))

## From https://www.bioinformatics.org/sms2/dna_mw.html (singled stranded linear).
## Units in Da.
synthgene16S_mol_wght_da<-50293.50
synthgeneITS_mol_wght_da<-51656.44

# Calculate number of synthgene molecules per sample.
## 16S. 1077665 molecules per sample.
(num_synthgenge_molecules_per_sample_16S<-0.09*(30110868216750/50)*(1/synthgene16S_mol_wght_da))
## ITS. 1049231 molecules per sample.
(num_synthgenge_molecules_per_sample_ITS<-0.09*(30110868216750/50)*(1/synthgeneITS_mol_wght_da))

# Get absolute read counts for sample types.
read_quantity_abs<-rowSums(tb16S_bio_scl)
read_quantity_abs<-data.frame(SampleID=names(read_quantity_abs),Quantity=read_quantity_abs)

# Add types to salamander read quantity.
read_quantity_abs$Type<-meta$Type

# Quick plot of absolute abundances.
library(ggplot2)
(g<-ggplot(data=read_quantity_abs,aes(x=Type,y=Quantity))+
    geom_violin(bw=0.375,trim=F)+
    geom_point(pch=21,alpha=0.5)+
    ylab("Number of DNA Molecules in Sample\nRelative to Number of Synthgene Molecules")+
    xlab("Sample Type")+
    ggtitle("Absolute Abundance by Sample Type",
            subtitle="16S Region")+
    theme_light()+
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          plot.subtitle=element_text(hjust=0.5),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
    scale_y_log10())

# Save absolute abundance plot.
ggsave(filename="Output_16S_2/DNARelative_SampleType_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# Convert scaled abundances to molecules of DNA per 15 uL DNA Extract.
tb16S_bio_scl<-tb16S_bio_scl*num_synthgenge_molecules_per_sample_16S

# Convert scaled abundances to molecules of DNA in sample.
# (15 uL of DNA extract out of 50 was used in each of the 16S and ITS processes.)
# (Upscale number of molecules to represent the number present in each sample.)
# (Units: swab, 500 mL water, 250 mg soil.)
tb16S_bio_scl<-tb16S_bio_scl*(50/15)

# Get absolute read counts for sample types.
read_quantity_abs<-rowSums(tb16S_bio_scl)
read_quantity_abs<-data.frame(SampleID=names(read_quantity_abs),Quantity=read_quantity_abs)

# Add types to salamander read quantity.
read_quantity_abs$Type<-meta$Type

# Quick plot of absolute abundances.
library(ggplot2)
(g<-ggplot(data=read_quantity_abs,aes(x=Type,y=Quantity))+
    geom_violin(bw=0.375,trim=F)+
    geom_point(pch=21,alpha=0.5)+
    ylab("Number of DNA Molecules in Sample")+
    xlab("Sample Type")+
    ggtitle("Absolute Abundance by Sample Type",
            subtitle="16S Region")+
    theme_light()+
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          plot.subtitle=element_text(hjust=0.5),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
    scale_y_log10())

# Save absolute abundance plot.
ggsave(filename="Output_16S_2/DNACount_SampleType_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# Check the positive mock community controls.

# Isolate the mock community controls.
mock<-tb16S_bio_props[row.names(tb16S_bio_props) %in% c("MC_1","MC_2"),]

# Remove unobserved taxa from the mock community controls.
mock<-mock[,colSums(mock)!=0]

# Get reformatted mock community taxa names.
split_names<-strsplit(colnames(mock),split="\\.")
colnames(mock)<-sapply(split_names,FUN=function(x){
    first_letters<-substr(x[-length(x)],start=1,stop=1)
    last_uppercase_level<-max(which(first_letters==toupper(first_letters)))
    paste(paste(x[last_uppercase_level:(length(x)-1)],collapse=" "),sapply(strsplit(x[length(x)],split="_"),"[[",2))
  }
)

# Get mock community sample names as a field.
mock$Sample<-row.names(mock)
row.names(mock)<-1:nrow(mock)

# Reshape mock community composition data.
mock<-reshape(data=mock,
        varying=colnames(mock)[-ncol(mock)],
        v.names="Proportional Abundance",
        timevar="Taxon",
        times=colnames(mock)[-ncol(mock)],
        new.row.names=1:prod(dim(mock[,-ncol(mock)])),
        direction="long")

# Remove the ID field.
mock<-mock[,-which(colnames(mock)=="id")]

# Quick plot of mock communities.
library(ggplot2)
(g<-ggplot(data=mock,aes(x=Sample,y=`Proportional Abundance`,fill=Taxon))+
    geom_bar(position="stack",stat="identity",color="black")+
    ylab("Proportional Abundance")+
    ggtitle("Proportional Abudances in Mock Community Samples",
            subtitle="16S Region")+
    theme_light()+
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          plot.subtitle=element_text(hjust=0.5),
          legend.title.align=0.5))

# Save mock community plot.
ggsave(filename="Output_16S_2/MockCommunitySamples_16S.tiff",plot=g,width=7,height=5,units="in",dpi=300)

# ZymoBIOMICS Microbial Community DNA Standard. Theoretical Composition Based on Genomic DNA: Listeria monocytogenes - 12%, Pseudomonas aeruginosa - 12%, Bacillus subtilis - 12%, Escherichia coli - 12%, Salmonella enterica - 12%, Lactobacillus fermentum - 12%, Enterococcus faecalis - 12%, Staphylococcus aureus - 12%, Saccharomyces cerevisiae - 2%, and Cryptococcus neoformans - 2%.

# From: https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-dna-standard

# Write out the checkThese table.
write.csv(checkThese,file="Output_16S_2/CheckTheseSamples_16S.csv",row.names=F)

# See Check_ASV_Table_16S.R for closing thoughts.

# Done!
print("Done!")
