setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs/FUNGuild")

FunGuild<- read.table("ASVTable.taxa.guilds.txt",sep="\t",header=TRUE)
length(unique(FunGuild$guild))
length(grep("Animal Pathogen",FunGuild$guild))
length(grep("Animal Endosymbiont",FunGuild$guild))
length(grep("na", FunGuild$guild))
#Of the 2099 ASVs, only 818 have any guild information, 1281 have no information
View(as.matrix(unique(FunGuild$guild)))

grep("Animal Pathogen", FunGuild$guild)
length(FunGuild$guild)

#I should make a key to know which fungi are animal pathogens and which are plant associated
FungusKey <- matrix(nrow=2099, ncol=2)
colnames(FungusKey)<- c("Animal Pathogen","Plant Associated")
FungusKey

FungusKey[grep("Animal Pathogen",FunGuild$guild),1]<- 1
FungusKey[grep("Plant Pathogen",FunGuild$guild),2]<- 1
FungusKey[grep("Endophyte",FunGuild$guild),2]<- 1
FungusKey[grep("Epiphyte",FunGuild$guild),2]<- 1
FungusKey[grep("Bryophyte Parasite",FunGuild$guild),2]<- 1
FungusKey[grep("Ectomycorrhizal",FunGuild$guild),2]<- 1
FungusKey[grep("Plant Parasite",FunGuild$guild),2]<- 1
FungusKey[is.na(FungusKey)]<-0
colSums(FungusKey)
length(which(rowSums(FungusKey)==2))
#There are 278 animal pathogens, and 216 of them are also plant associated in some way
#There are 456 plant associated fungi, and 216 of them are also animal pathogens
#That means that these two categories cover 518 of the 818 ASVs with known guilds



#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Loading in relevant packages
library(qiime2R)#We need this to bring in the qza objects
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(stringr)#For strings
library(lme4)#For mixed effect models
library(car)
library(agricolae)

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

#I looked at the rarefaction curve in Qiime2 and want to hit 5,000 reads
#I have a number of samples with double digits or less, so they go regardless
#I need to remove things with too few reads myself, lose 7 samples
MetadataReduced <- metadata[-which(colSums(FungTable)<5000),]
FungTableReduced <- FungTable[,-which(colSums(FungTable)<5000)]
dim(FungTableReduced)

#I need to make a table of just the inclusive plant associated fungi
PLantAssocTable <- FungTableReduced[which(FungusKey[,2]==1),]

#I need an exclusive animal pathogen table
AniPathTable <- FungTableReduced[which(FungusKey[,1] == 1 & FungusKey[,2] == 0),]

#Pulling in diet stuff
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

sum(table(Plant$SCIENTIFIC_NAME))
sum(is.na(Plant$SCIENTIFIC_NAME))
#A vast majority of these (5,313) have some level of identification
#Only nine of them are completely unidentified

#Let's just get these 60 samples from the plants
PlantsReduced <- Plant[,which(MetadataReduced$sampleid[1]==colnames(Plant))]
SelectListPlant <- which(MetadataReduced$sampleid[1]==colnames(Plant))
for(i in 2:length(MetadataReduced$sampleid)){
  PlantsReduced <- cbind(PlantsReduced,Plant[,which(MetadataReduced$sampleid[i]==colnames(Plant))])
  SelectListPlant <- c(SelectListPlant,which(MetadataReduced$sampleid[i]==colnames(Plant)))
}

dim(PlantsReduced)

colnames(PlantsReduced) <- colnames(Plant)[SelectListPlant]


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



#Now for mantel tests
library(ecodist)
#Let's make a subsample for mantel tests
#Now I want to create a single subsample for beta diversity analyses
set.seed(1210)

Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
dim(Subsample)
rowSums(Subsample)
#I should create an Animal Pathogen subsample
APSubsample <- Subsample[,which(FungusKey[,1]==1 & FungusKey[,2]==0)]
dim(APSubsample)
#Now a plant associated subsample
PASubsample <- Subsample[,which(FungusKey[,2]==1)]
dim(PASubsample)


ecodist::mantel(formula=(vegdist((Subsample), method="bray", binary=FALSE)) ~(vegdist(t(InvertsReduced), method="bray", binary=FALSE)), nperm=100000)

ecodist::mantel(formula=(vegdist((Subsample), method="bray", binary=FALSE)) ~(vegdist(t(PlantsReduced), method="bray", binary=FALSE)), nperm=100000)

#Some samples don't have plant associated taxa, so those samples must go
#The mantel needs the same samples in both matrices
PASubSample2 <- PASubsample[-which(rowSums(PASubsample)==0),]
InvertsReduced2 <- InvertsReduced[,-which(rowSums(PASubsample)==0)]
PlantsReduced2 <- PlantsReduced[,-which(rowSums(PASubsample)==0)]

ecodist::mantel(formula=(vegdist((PASubSample2), method="bray", binary=FALSE)) ~(vegdist(t(InvertsReduced2), method="bray", binary=FALSE)), nperm=100000)

ecodist::mantel(formula=(vegdist((PASubSample2), method="bray", binary=FALSE)) ~(vegdist(t(PlantsReduced2), method="bray", binary=FALSE)), nperm=100000)

#Let's drop samples that lack animal pathogens
APSubsample2 <- APSubsample[-which(rowSums(APSubsample)==0),]
InvertsReduced3 <- InvertsReduced[,-which(rowSums(APSubsample)==0)]
PlantsReduced3 <- PlantsReduced[,-which(rowSums(APSubsample)==0)]

ecodist::mantel(formula=(vegdist((APSubsample2), method="bray", binary=FALSE)) ~(vegdist(t(InvertsReduced3), method="bray", binary=FALSE)), nperm=100000)

ecodist::mantel(formula=(vegdist((APSubsample2), method="bray", binary=FALSE)) ~(vegdist(t(PlantsReduced3), method="bray", binary=FALSE)), nperm=100000)
