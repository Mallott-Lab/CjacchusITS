#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData")

#Loading relevant packages
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(qiime2R)#We need this to bring in the qza objects
library(stringr)#For strings

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

#It doesn't look like they rarefied the diet reads in the other manuscript, so I won't either for consistency
#I am going to to calculate beta diversity for the two diets

#I need my metadata
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

#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Let's bring in the objects from Qiime2
tabletrimmedQiime2Object<-read_qza("table-trimmed.qza", rm = FALSE)
taxtableQiime2Object<-read_qza("taxonomy.qza")

#Let's pull out what I specifically want to work with
FungTable<-tabletrimmedQiime2Object$data
TaxTable<-taxtableQiime2Object$data
#I looked at the rarefaction curve in Qiime2 and want to hit 5,000 reads
#I have a number of samples with double digits or less, so they go regardless
#I need to remove things with too few reads myself, lose 7 samples
MetadataReduced <- metadata[-which(colSums(FungTable)<5000),]

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

#There are approximately triple the number of samples that the mycobiome has at this point
#I am going to pull out the plant information for the samples that I have in my mycobiome work

PlantsReduced <- Plant[,which(MetadataReduced$sampleid[1]==colnames(Plant))]
SelectListPlant <- which(MetadataReduced$sampleid[1]==colnames(Plant))
for(i in 2:length(MetadataReduced$sampleid)){
  PlantsReduced <- cbind(PlantsReduced,Plant[,which(MetadataReduced$sampleid[i]==colnames(Plant))])
  SelectListPlant <- c(SelectListPlant,which(MetadataReduced$sampleid[i]==colnames(Plant)))
}

dim(PlantsReduced)

colnames(PlantsReduced) <- colnames(Plant)[SelectListPlant]

#There are some taxa with no reads
#I want to cut them out
#I am also getting taxonomy ready
PlantTaxonomy<-Plant[,1:2]
PlantTaxonomy<-PlantTaxonomy[-which(rowSums(PlantsReduced)==0),]
PlantsReduced<-PlantsReduced[-which(rowSums(PlantsReduced)==0),]

#Let's get my set of samples for the invertebrates
InvertsReduced <- Invert[,which(MetadataReduced$sampleid[1]==colnames(Invert))]
SelectListInvert <- which(MetadataReduced$sampleid[1]==colnames(Invert))
for(i in 2:length(MetadataReduced$sampleid)){
  InvertsReduced <- cbind(InvertsReduced,Invert[,which(MetadataReduced$sampleid[i]==colnames(Invert))])
  SelectListInvert <- c(SelectListInvert,which(MetadataReduced$sampleid[i]==colnames(Invert)))
}

dim(InvertsReduced)
colnames(InvertsReduced) <- colnames(Invert)[SelectListInvert]

#I have some empty rows in the invert table
InvertTaxonomy<-Invert[,1:2]
InvertTaxonomy<-InvertTaxonomy[-which(rowSums(InvertsReduced)==0),]
InvertsReduced<-InvertsReduced[-which(rowSums(InvertsReduced)==0),]
#There are 410 inverts left


#I am going to create a Bray-Curtis dissimilarity matrix
InvertBrayDist<-vegdist(t(InvertsReduced),method="bray",binary=FALSE)
#Here is a Jaccard distance matrix
InvertJaccDist<- vegdist(t(InvertsReduced),method="jaccard",binary=TRUE)

#Now the data are ready for analyses
#Let's start with the PERMANOVA of the Bray-Curtis dissimilarity matrix
#Dropped Time from the PERMANOVA because it is nested in season and PERMANOVA doesn't quite handle that well
InvertBrayPERMANOVA<-adonis2(InvertBrayDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

InvertBrayPERMANOVA
#Season is significant in the PERMANOVA

#Let's run NMDS
InvertBrayNMDS<-metaMDS(InvertBrayDist, k = 2,autotransform=FALSE,try=100,trymax=100)
InvertBrayNMDS

#Let's plot season it is significant
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(InvertBrayNMDS$points,pch=19, cex= 0.5,col="white", main= "Bray-Curtis NMDS")

#Seasons were significant, so let's mark those with colors
points(InvertBrayNMDS$points[which(MetadataReduced$Season=="Dry"),], pch = 19, cex = 1.5, col = "#E69F00")
points(InvertBrayNMDS$points[which(MetadataReduced$Season=="Wet"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(InvertBrayNMDS$points,MetadataReduced$Season,col=c("#E69F00","#56B4E9"))
#There is one crazy sample there

#Now to look at Jaccard
InvertJaccardPERMANOVA<-adonis2(InvertJaccDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

InvertJaccardPERMANOVA
#Season is significant

#Let's run NMDS
InvertJaccardNMDS<-metaMDS(InvertJaccDist, k = 2,autotransform=FALSE,try=100,trymax=100)
InvertJaccardNMDS

#Let's plot season it is significant
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(InvertJaccardNMDS$points,pch=19, cex= 0.5,col="white", main= "Jaccard NMDS")

#Seasons were significant, so let's mark those with colors
points(InvertJaccardNMDS$points[which(MetadataReduced$Season=="Dry"),], pch = 19, cex = 1.5, col = "#E69F00")
points(InvertJaccardNMDS$points[which(MetadataReduced$Season=="Wet"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(InvertJaccardNMDS$points,MetadataReduced$Season,col=c("#E69F00","#56B4E9"))
#There is one crazy sample there


#I am going to create a Bray-Curtis dissimilarity matrix
PlantBrayDist<-vegdist(t(PlantsReduced),method="bray",binary=FALSE)
#Here is a Jaccard distance matrix
PlantJaccDist<- vegdist(t(PlantsReduced),method="jaccard",binary=TRUE)

#Now the data are ready for analyses
#Let's start with the PERMANOVA of the Bray-Curtis dissimilarity matrix
#Dropped Time from the PERMANOVA because it is nested in season and PERMANOVA doesn't quite handle that well
PlantBrayPERMANOVA<-adonis2(PlantBrayDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

PlantBrayPERMANOVA
#Preservative is significant in the PERMANOVA

#Let's run NMDS
PlantBrayNMDS<-metaMDS(PlantBrayDist, k = 2,autotransform=FALSE,try=100,trymax=100)
PlantBrayNMDS

#Let's plot season it is significant
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(PlantBrayNMDS$points,pch=19, cex= 0.5,col="white", main= "Bray-Curtis NMDS")

#Preservative was significant, so let's mark those with colors
points(PlantBrayNMDS$points[which(MetadataReduced$Preservative=="RNAlater"),], pch = 19, cex = 1.5, col = "#E69F00")
points(PlantBrayNMDS$points[which(MetadataReduced$Preservative=="Ethanol"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(PlantBrayNMDS$points,MetadataReduced$Preservative,col=c("#E69F00","#56B4E9"))
#There is one crazy sample there

#Now to look at Jaccard
PlantJaccardPERMANOVA<-adonis2(PlantJaccDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

PlantJaccardPERMANOVA
#Preservative is significant

#Let's run NMDS
PlantJaccardNMDS<-metaMDS(PlantJaccDist, k = 2,autotransform=FALSE,try=100,trymax=100)
PlantJaccardNMDS

#Let's plot season it is significant
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(PlantJaccardNMDS$points,pch=19, cex= 0.5,col="white", main= "Jaccard NMDS")

#Preservative was significant, so let's mark those with colors
points(PlantJaccardNMDS$points[which(MetadataReduced$Preservative=="RNAlater"),], pch = 19, cex = 1.5, col = "#E69F00")
points(PlantJaccardNMDS$points[which(MetadataReduced$Preservative=="Ethanol"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(PlantJaccardNMDS$points,MetadataReduced$Preservative,col=c("#E69F00","#56B4E9"))
#There is one crazy sample there