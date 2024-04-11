library(ccrepe)

#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Loading in relevant packages
library(qiime2R)#We need this to bring in the qza objects
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(stringr)#For strings
library(lme4)#For mixed effect models

#Let's bring in the objects from Qiime2
tabletrimmedQiime2Object<-read_qza("table-trimmed.qza", rm = FALSE)
taxtableQiime2Object<-read_qza("taxonomy.qza")

#Let's pull out what I specifically want to work with
FungTable<-tabletrimmedQiime2Object$data
TaxTable<-taxtableQiime2Object$data
#I put these files out to run things through FUNGuild
#write.table(TaxTable,"TaxTable.txt", sep = "\t")
#I had to delete the column of index numbers in this
#write.table(FungTable,"ASVTable.txt", sep = "\t")

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
#Let's see how much the read depth varies
barplot(colSums(FungTable))
barplot(colSums(FungTable),ylim=c(0,10000))

#I looked at the rarefaction curve in Qiime2 and want to hit 5,000 reads
#I have a number of samples with double digits or less, so they go regardless
#I need to remove things with too few reads myself, lose 7 samples
MetadataReduced <- metadata[-which(colSums(FungTable)<5000),]
FungTableReduced <- FungTable[,-which(colSums(FungTable)<5000)]
dim(FungTableReduced)

#I need to remove samples that don't have sex information
dim(FungTableReduced)
FungTableReduced <- FungTableReduced[,-which(MetadataReduced$Sex=="Unknown")]
dim(FungTableReduced)

#I should drop them from the metadata
MetadataReduced <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced)

#I should make some of these into factors
MetadataReduced$Sex <-as.factor(MetadataReduced$Sex)
MetadataReduced$Group <-as.factor(MetadataReduced$Group)
MetadataReduced$Season <-as.factor(MetadataReduced$Season)
MetadataReduced$Reproductive <-as.factor(MetadataReduced$Reproductive)
MetadataReduced$Individual <-as.factor(MetadataReduced$Individual)
MetadataReduced$Preservative <-as.factor(MetadataReduced$Preservative)
MetadataReduced$Age <-as.factor(MetadataReduced$Age)
MetadataReduced$Time2 <-as.factor(MetadataReduced$Time2)

#Now I want to create a single subsample for beta diversity analyses
set.seed(1210)

Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
dim(Subsample)
rowSums(Subsample)
#I should remove any ASVs that disappeared as I removed samples and subsampled
dim(Subsample)
Subsample <- Subsample[,-which(colSums(Subsample)==0)]
dim(Subsample)
Subsample



#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData")

#Let's load in the necessary packages
library(vegan)#distances and diversity metrics
library(fossil)#For Chao1 diversity
#Now I am loading the in the data

#This is the invert_results_clean to see which invertebrate taxa are detected in the marmoset samples
Invert <- read.table("invert_results_clean.txt",sep = "\t", header=TRUE)

sum(table(Invert$TaxID))
sum(is.na(Invert$TaxID))
#It looks like 972 of the taxa are unidentified in any way and 396 have some identification

#Let's load in the plant taxa to see which plants are detected in the marmoset samples
Plant <- read.table("plant_results_clean.txt",sep = "\t", header=TRUE)




#Let's get my set of samples for the invertebrates
InvertsReduced <- Invert[,which(MetadataReduced$sampleid[1]==colnames(Invert))]
SelectListInvert <- which(MetadataReduced$sampleid[1]==colnames(Invert))
for(i in 2:length(MetadataReduced$sampleid)){
  InvertsReduced <- cbind(InvertsReduced,Invert[,which(MetadataReduced$sampleid[i]==colnames(Invert))])
  SelectListInvert <- c(SelectListInvert,which(MetadataReduced$sampleid[i]==colnames(Invert)))
}

dim(InvertsReduced)
colnames(InvertsReduced) <- colnames(Invert)[SelectListInvert]
InvertsReduced



PlantsReduced <- Plant[,which(MetadataReduced$sampleid[1]==colnames(Plant))]
SelectListPlant <- which(MetadataReduced$sampleid[1]==colnames(Plant))
for(i in 2:length(MetadataReduced$sampleid)){
  PlantsReduced <- cbind(PlantsReduced,Plant[,which(MetadataReduced$sampleid[i]==colnames(Plant))])
  SelectListPlant <- c(SelectListPlant,which(MetadataReduced$sampleid[i]==colnames(Plant)))
}

dim(PlantsReduced)

colnames(PlantsReduced) <- colnames(Plant)[SelectListPlant]

#I need relative abundance of Fungi, Plants, and Invertebrates for ccrepe
RelPlants <- matrix(nrow=dim(PlantsReduced)[1],ncol=dim(PlantsReduced)[2])
for(i in 1:dim(PlantsReduced)[1]){
  for(j in 1:dim(PlantsReduced)[2]){
    RelPlants[i,j] <- PlantsReduced[i,j]/colSums(PlantsReduced)[j]
  }
}
colSums(RelPlants)


#I need relative abundance of Fungi, Plants, and Invertebrates for ccrepe
RelInverts <- matrix(nrow=dim(InvertsReduced)[1],ncol=dim(InvertsReduced)[2])
for(i in 1:dim(InvertsReduced)[1]){
  for(j in 1:dim(InvertsReduced)[2]){
    RelInverts[i,j] <- InvertsReduced[i,j]/colSums(InvertsReduced)[j]
  }
}
colSums(RelInverts)
dim(RelInverts)
RelInverts<-RelInverts[-which(rowSums(RelInverts)==0),]
dim(RelInverts)

RelFung <- matrix(nrow=dim(FungTableReduced)[1],ncol=dim(FungTableReduced)[2])
for(i in 1:dim(FungTableReduced)[1]){
  for(j in 1:dim(FungTableReduced)[2]){
    RelFung[i,j] <- FungTableReduced[i,j]/colSums(FungTableReduced)[j]
  }
}
colSums(RelFung)
dim(RelFung)
RelFung<-RelFung[-which(rowSums(RelFung)==0),]
dim(RelFung)

colSums(RelPlants)
dim(RelPlants)
RelPlants<-RelPlants[-which(rowSums(RelPlants)==0),]
dim(RelPlants)

#I am going to set a random seed for consistency
set.seed(121023)
FungInvCCREPE <- ccrepe(x = t(RelFung), y = t(RelInverts),verbose=TRUE, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"),min.subj = 26, make.output.table=TRUE)
#minimum subjects should be at least 50% from Sharma et al 2022
#I started this at 13:22 and finished around
#There were 50 or more warnings, all saying that the standard deviation is zero
#It did work though
#The Sharma paper highlighted sim scores with absolute values above 0.6 and q values < 0.01

#I want to put the results into documents so I don't have to run this again
#setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/CCREPE_RESULTS")
#FungInvCCREPE_pvalues <- FungInvCCREPE$p.values
#FungInvCCREPE_qvalues <- FungInvCCREPE$q.values
#FungInvCCREPE_zstat <- FungInvCCREPE$z.stat
#FungInvCCREPE_simscore <- FungInvCCREPE$sim.score
#FungInvCCREPE_table <- FungInvCCREPE$output.table
#write.table(FungInvCCREPE_pvalues, file = "FungInvCCREPE_pvalues.txt", sep = "\t")
#write.table(FungInvCCREPE_qvalues, file = "FungInvCCREPE_qvalues.txt", sep = "\t")
#write.table(FungInvCCREPE_zstat, file = "FungInvCCREPE_zstat.txt", sep = "\t")
#write.table(FungInvCCREPE_simscore, file = "FungInvCCREPE_simscore.txt", sep = "\t")
#write.table(FungInvCCREPE_table, file = "FungInvCCREPE_table.txt", sep = "\t")

#I need to format these differently for Cytoscape
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/CCREPE_RESULTS")

FungInvCCREPE_pvalues <- read.table(file = "FungInvCCREPE_pvalues.txt", sep = "\t")
FungInvCCREPE_qvalues <- read.table(file = "FungInvCCREPE_qvalues.txt", sep = "\t")
FungInvCCREPE_zstat <- read.table(file = "FungInvCCREPE_zstat.txt", sep = "\t")
FungInvCCREPE_simscore <- read.table(file = "FungInvCCREPE_simscore.txt", sep = "\t")

dim(FungInvCCREPE_pvalues)
dim(RelInverts)#I have 410 invertebrates
dim(RelFung)#I have 1839 fungi
#I need to make a new thing

#Most of the taxa don't interact at all
#Of the 753990 potential interactions, there are 96
CytoscapeData <- matrix(nrow = 96, ncol = 6)


CytoscapeData[,3]<-as.vector(as.matrix(FungInvCCREPE_pvalues))[-which(is.na(as.vector(as.matrix(FungInvCCREPE_pvalues))))]
CytoscapeData[,4]<-as.vector(as.matrix(FungInvCCREPE_qvalues))[-which(is.na(as.vector(as.matrix(FungInvCCREPE_qvalues))))]
CytoscapeData[,5]<-as.vector(as.matrix(FungInvCCREPE_simscore))[-which(is.na(as.vector(as.matrix(FungInvCCREPE_simscore))))]
CytoscapeData[,6]<-as.vector(as.matrix(FungInvCCREPE_zstat))[-which(is.na(as.vector(as.matrix(FungInvCCREPE_zstat))))]
InvertList <- Invert$ID
InvertList <- InvertList[-which(rowSums(InvertsReduced)==0)]
FungTableReduced <- FungTableReduced[-which(rowSums(FungTableReduced)==0),]
#Now I just need the fungus and invertebrate labels
COUNT<-0
for(i in 1:dim(RelInverts)[1]){
  for(j in 1:dim(RelFung)[1]){
    if(is.na(FungInvCCREPE_pvalues[j,i])==FALSE){COUNT<-COUNT+1}
    if(is.na(FungInvCCREPE_pvalues[j,i])==FALSE){CytoscapeData[COUNT,1]<-rownames(FungTableReduced)[j]}
    if(is.na(FungInvCCREPE_pvalues[j,i])==FALSE){CytoscapeData[COUNT,2]<-InvertList[i]}
  }
}
#The fungus labels worked
#write.table(CytoscapeData, file = "CytoscapeData.txt", sep = "\t")



#This function frees up memory without removing objects
#gc()
#set.seed(11623)
#FungPlantCCREPE <- ccrepe(x = t(RelFung), y = t(RelPlants),verbose=TRUE, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"),min.subj = 26)
#started at 12:43, failed at 13:02

#setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/CCREPE_RESULTS")
#FungPlantCCREPE_pvalues <- FungPlantCCREPE$p.values
#FungPlantCCREPE_qvalues <- FungPlantCCREPE$q.values
#FungPlantCCREPE_zstat <- FungPlantCCREPE$z.stat
#FungPlantCCREPE_simscore <- FungPlantCCREPE$sim.score
#write.table(FungPlantCCREPE_pvalues, file = "FungPlantCCREPE_pvalues.txt", sep = "\t")
#write.table(FungPlantCCREPE_qvalues, file = "FungPlantCCREPE_qvalues.txt", sep = "\t")
#write.table(FungPlantCCREPE_zstat, file = "FungPlantCCREPE_zstat.txt", sep = "\t")
#write.table(FungPlantCCREPE_simscore, file = "FungPlantCCREPE_simscore.txt", sep = "\t")



#I am going to read in the results of the CCREPE and format them for Cytoscape

setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/CCREPE_RESULTS")

#FungPlantCCREPE_pvalues <- read.table(file = "FungPlantCCREPE_pvalues.txt", sep = "\t")
#FungPlantCCREPE_qvalues <- read.table(file = "FungPlantCCREPE_qvalues.txt", sep = "\t")
#FungPlantCCREPE_zstat <- read.table(file = "FungPlantCCREPE_zstat.txt", sep = "\t")
#FungPlantCCREPE_simscore <- read.table(file = "FungPlantCCREPE_simscore.txt", sep = "\t")
#CytoscapePlantData <- read.table(file = "CytoscapeDataPlants.txt", sep = "\t",header=TRUE)
#CytoscapeInvertData <- read.table(file = "CytoscapeDataInverts.txt", sep = "\t",header=TRUE)


#Most of the taxa don't interact at all
#Of the 3185148 potential interactions, there are 3185004 NAs
#There are only 144 that are not NA
CytoscapeData2 <- matrix(nrow = 144, ncol = 6)


CytoscapeData2[,3]<-as.vector(as.matrix(FungPlantCCREPE_pvalues))[-which(is.na(as.vector(as.matrix(FungPlantCCREPE_pvalues))))]
CytoscapeData2[,4]<-as.vector(as.matrix(FungPlantCCREPE_qvalues))[-which(is.na(as.vector(as.matrix(FungPlantCCREPE_qvalues))))]
CytoscapeData2[,5]<-as.vector(as.matrix(FungPlantCCREPE_simscore))[-which(is.na(as.vector(as.matrix(FungPlantCCREPE_simscore))))]
CytoscapeData2[,6]<-as.vector(as.matrix(FungPlantCCREPE_zstat))[-which(is.na(as.vector(as.matrix(FungPlantCCREPE_zstat))))]
PlantList <- Plant$ID
PlantList <- PlantList[-which(rowSums(PlantsReduced)==0)]
FungTableReduced <- FungTableReduced[-which(rowSums(FungTableReduced)==0),]
#Now I just need the fungus and invertebrate labels
COUNT<-0
for(i in 1:dim(RelPlants)[1]){
  for(j in 1:dim(RelFung)[1]){
    if(is.na(FungPlantCCREPE_pvalues[j,i])==FALSE){COUNT<-COUNT+1}
    if(is.na(FungPlantCCREPE_pvalues[j,i])==FALSE){CytoscapeData2[COUNT,1]<-rownames(FungTableReduced)[j]}
    if(is.na(FungPlantCCREPE_pvalues[j,i])==FALSE){CytoscapeData2[COUNT,2]<-PlantList[i]}
  }
}
#The fungus labels worked
#write.table(CytoscapeData2, file = "CytoscapeDataPlants.txt", sep = "\t")
