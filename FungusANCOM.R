library(ANCOMBC)
library(phyloseq)

#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Loading in relevant packages
library(qiime2R)#We need this to bring in the qza objects
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(stringr)#For strings
library(lme4)#For mixed effect models
library(car)#for type III ANOVA
#library(agricolae)#for the Tukey HSD

#Let's bring in the objects from Qiime2
tabletrimmedQiime2Object<-read_qza("table-trimmed.qza", rm = FALSE)
taxtableQiime2Object<-read_qza("taxonomy.qza")

#Let's pull out what I specifically want to work with
FungTable<-tabletrimmedQiime2Object$data
TaxTable<-taxtableQiime2Object$data

#Now I want to bring in my metadata
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork")
metadata<- read.table("Cjacchus_ITS_metadata.txt",sep="\t",header=TRUE)
#This metadata file is not in the same sequence as the table
metadata<- metadata[order(metadata$sampleid),]
#I want to make a change to reproductive status to distinguish non-pregnant females from males
metadata$Reproductive[which(metadata$Sex=="Male")]<-"NPMale"
table(metadata$Reproductive)
#I want to drop a bunch of the empty columns
dim(metadata)
metadata<- metadata[,-c(7,8,10,15,17,18,19,23)]
dim(metadata)
#I want to make a Season/year column
metadata$Time2<-str_c(metadata$Season,str_split_i(metadata$Date,"/",3))
#I also want to combine group and individual names
#Some of the individual names are reused between groups
metadata$Individual
metadata$Individual[which(metadata$Individual!="")]<-str_c(metadata$Individual[which(metadata$Individual!="")],metadata$Group[which(metadata$Individual!="")])
metadata$Individual

#I need to remove things with too few reads myself, lose 7 samples
MetadataReduced <- metadata[-which(colSums(FungTable)<5000),]
FungTableReduced <- FungTable[,-which(colSums(FungTable)<5000)]
dim(FungTableReduced)

FungTableReduced <- FungTableReduced[,-which(MetadataReduced$Sex=="Unknown")]
dim(FungTableReduced)
MetadataReduced <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced)

#Let's get rid of things that no longer exist
length(which(rowSums(FungTableReduced)==0))

TaxTableReduced <- TaxTable[-which(rowSums(FungTableReduced)==0),]
dim(TaxTableReduced)
rownames(TaxTableReduced) <- TaxTableReduced[,1]
FungTableReduced <- FungTableReduced[-which(rowSums(FungTableReduced)==0),]
dim(FungTableReduced)

rownames(MetadataReduced) <- MetadataReduced$sampleid

#Now let's make a phyloseq object
FungPhylo <- phyloseq(otu_table(FungTableReduced,taxa_are_rows = TRUE),tax_table(as.matrix(TaxTableReduced)),sample_data(MetadataReduced))

#I need to make a tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(FungPhylo)

FungANCOM <- ancombc2(data = tse, fix_formula = "Sex + Season", p_adj_method = "BH")
FungANCOM$res
View(FungANCOM$res)
FungANCOM$res[which(FungANCOM$res$diff_SexMale==TRUE),]
FungANCOM$res[which(FungANCOM$res$diff_SeasonWet==TRUE),]
FungANCOM$res$lfc_SexMale

SexNames<-c("Fungus_2", "Nothophoma", "Cladosporium_1", "Leptosillia","Cladosporium_2","Fungus_8","Fungus_9")

par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
SexLFC <- barplot(FungANCOM$res$lfc_SexMale[which(FungANCOM$res$diff_SexMale==TRUE)],ylim=c(-3.5,3),names=SexNames, ylab = "Logfold Change",xlab = "Taxon", main = "Differences by Sex", col = c("red","blue","green","orange","hotpink", "gray","brown"))

arrows(x0=SexLFC, y0 = FungANCOM$res$lfc_SexMale[which(FungANCOM$res$diff_SexMale==TRUE)] + FungANCOM$res$se_SexMale[which(FungANCOM$res$diff_SexMale==TRUE)], y1 =FungANCOM$res$lfc_SexMale[which(FungANCOM$res$diff_SexMale==TRUE)] - FungANCOM$res$se_SexMale[which(FungANCOM$res$diff_SexMale==TRUE)], angle = 90, code = 3, length = 0.1)

SeasonNames<-c("Nothophoma","Cladosporium_1","Cladosporium_2","Candida parapsilosis")

SeasonLFC <- barplot(FungANCOM$res$lfc_SeasonWet[which(FungANCOM$res$diff_SeasonWet==TRUE)],ylim=c(-6,5), ylab = "Logfold Change",xlab = "Taxon", names = SeasonNames,main = "Differences by Season", col = c("blue","green","hotpink","purple"))

arrows(x0=SeasonLFC, y0 = FungANCOM$res$lfc_SeasonWet[which(FungANCOM$res$diff_SeasonWet==TRUE)] + FungANCOM$res$se_SeasonWet[which(FungANCOM$res$diff_SeasonWet==TRUE)], y1 =FungANCOM$res$lfc_SeasonWet[which(FungANCOM$res$diff_SeasonWet==TRUE)] - FungANCOM$res$se_SeasonWet[which(FungANCOM$res$diff_SeasonWet==TRUE)], angle = 90, code = 3, length = 0.1)
