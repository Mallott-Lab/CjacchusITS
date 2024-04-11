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

#It doesn't look like they rarefied the diet reads in the other manuscript, so I won't either for consistency
#I am going to calculate alpha diversity metrics for the diet (Shannon, Simpson, Chao1, Pielou, and Richness)
#I am going to do this for each section of the diet (plants and invertebrates)

#I want to see how these diet based alpha diversity metrics interact with the fungal alpha diversity metrics
#I need to pull in the mycobiome data and calculate those metrics
#I ran lines 1 through 88 of the AlphaDiversityMixedEffectLinearModelsMarmosetMycobiome.R code file
#I am going to bring myself back to the 52 samples that are in the other analyses
dim(RarefiedAlphaDiversityTable)
RarefiedAlphaDiversityTable <- RarefiedAlphaDiversityTable[,-which(MetadataReduced$Sex == "Unknown")]
dim(RarefiedAlphaDiversityTable)

dim(MetadataReduced)
MetadataReduced <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced)

#I need to drop samples from my mycobiome alpha diversity metrics too

#I need to match these up to the diet data
MetadataReduced$sampleid
colnames(Plant)[3:180]
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

#I am going to make the alpha diversity metrics for these samples
#Here are alpha diversity metrics for the plants
ShannonDiversityPlants <- diversity(t(PlantsReduced), index="shannon")
SimpsonDiversityPlants <- diversity(t(PlantsReduced), index="simpson")
Chao1DiversityPlants <- apply(t(PlantsReduced), MARGIN = 1,chao1, taxa.row=FALSE)
RichnessPlants <- rowSums(ifelse(t(PlantsReduced)>0,1,0))
PielouPlants <- ShannonDiversityPlants/log(RichnessPlants)

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

#I am going to make the alpha diversity metrics for these samples
#Here are alpha diversity metrics for the invertebrates
ShannonDiversityInverts <- diversity(t(InvertsReduced), index="shannon")
SimpsonDiversityInverts <- diversity(t(InvertsReduced), index="simpson")
Chao1DiversityInverts <- apply(t(InvertsReduced), MARGIN = 1,chao1, taxa.row=FALSE)
RichnessInverts <- rowSums(ifelse(t(InvertsReduced)>0,1,0))
PielouInverts <- ShannonDiversityInverts/log(RichnessInverts)

#Shannon 
plot(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTable[1,], main = "Shannon", xlab = "Invertebrate Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTable[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTable[1,], main = "Shannon", xlab = "Plant Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTable[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Richness 
plot(x = RichnessInverts, y = RarefiedAlphaDiversityTable[2,], main = "Richness", xlab = "Invertebrate Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = RichnessInverts, y = RarefiedAlphaDiversityTable[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessPlants, y = RarefiedAlphaDiversityTable[2,], main = "Richness", xlab = "Plant Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = RichnessPlants, y = RarefiedAlphaDiversityTable[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Chao1 
plot(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTable[3,], main = "Chao1", xlab = "Invertebrate Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTable[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTable[3,], main = "Chao1", xlab = "Plant Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTable[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Simpson 
plot(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTable[4,], main = "Simpson", xlab = "Invertebrate Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTable[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTable[4,], main = "Simpson", xlab = "Plant Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTable[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Pielou 
plot(x = PielouInverts, y = RarefiedAlphaDiversityTable[5,], main = "Pielou's Evenness", xlab = "Invertebrate Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = PielouInverts, y = RarefiedAlphaDiversityTable[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouPlants, y = RarefiedAlphaDiversityTable[5,], main = "Pielou's Evenness", xlab = "Plant Diversity", ylab = "Mycobiome Diversity", pch = 19)
cor.test(x = PielouPlants, y = RarefiedAlphaDiversityTable[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)




#If I run lines 1 to 193 of the FUNGUILD code I can get the RarefiedAlphaDiversityTableBothAsAniPath and RarefiedAlphaDiversityTableBothAsPlant
#This can be run with the correlation tests like the overal mycobiome alpha diversity metrics


#Shannon 
plot(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[1,], main = "Shannon", xlab = "Invertebrate Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[1,], main = "Shannon", xlab = "Plant Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[1,], main = "Shannon", xlab = "Invertebrate Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[1,], main = "Shannon", xlab = "Plant Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Richness 
plot(x = RichnessInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[2,], main = "Richness", xlab = "Invertebrate Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = RichnessInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[2,], main = "Richness", xlab = "Plant Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = RichnessPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessInverts, y = RarefiedAlphaDiversityTableBothAsPlant[2,], main = "Richness", xlab = "Invertebrate Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = RichnessInverts, y = RarefiedAlphaDiversityTableBothAsPlant[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessPlants, y = RarefiedAlphaDiversityTableBothAsPlant[2,], main = "Richness", xlab = "Plant Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = RichnessPlants, y = RarefiedAlphaDiversityTableBothAsPlant[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Chao1 
plot(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[3,], main = "Chao1", xlab = "Invertebrate Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[3,], main = "Chao1", xlab = "Plant Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[3,], main = "Chao1", xlab = "Invertebrate Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[3,], main = "Chao1", xlab = "Plant Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Simpson 
plot(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[4,], main = "Simpson", xlab = "Invertebrate Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[4,], main = "Simpson", xlab = "Plant Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[4,], main = "Simpson", xlab = "Invertebrate Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableBothAsPlant[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[4,], main = "Simpson", xlab = "Plant Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableBothAsPlant[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Pielou 
plot(x = PielouInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[5,], main = "Pielou's Evenness", xlab = "Invertebrate Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = PielouInverts, y = RarefiedAlphaDiversityTableBothAsAniPath[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[5,], main = "Pielou's Evenness", xlab = "Plant Diversity", ylab = "Animal Pathogen Diversity", pch = 19)
cor.test(x = PielouPlants, y = RarefiedAlphaDiversityTableBothAsAniPath[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouInverts, y = RarefiedAlphaDiversityTableBothAsPlant[5,], main = "Pielou's Evenness", xlab = "Invertebrate Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = PielouInverts, y = RarefiedAlphaDiversityTableBothAsPlant[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouPlants, y = RarefiedAlphaDiversityTableBothAsPlant[5,], main = "Pielou's Evenness", xlab = "Plant Diversity", ylab = "Plant Associated Diversity", pch = 19)
cor.test(x = PielouPlants, y = RarefiedAlphaDiversityTableBothAsPlant[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)


#To get the plant associated mycobiome alpha diversity values run lines 193 to 222 in FUNGUILDCODE
RarefiedAlphaDiversityTableOnlyPlant


#Shannon 
plot(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[1,], main = "Shannon", xlab = "Invertebrate Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[1,], main = "Shannon", xlab = "Plant Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Richness 
plot(x = RichnessInverts, y = RarefiedAlphaDiversityTableOnlyPlant[2,], main = "Richness", xlab = "Invertebrate Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = RichnessInverts, y = RarefiedAlphaDiversityTableOnlyPlant[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessPlants, y = RarefiedAlphaDiversityTableOnlyPlant[2,], main = "Richness", xlab = "Plant Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = RichnessPlants, y = RarefiedAlphaDiversityTableOnlyPlant[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Chao1 
plot(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[3,], main = "Chao1", xlab = "Invertebrate Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[3,], main = "Chao1", xlab = "Plant Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Simpson 
plot(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[4,], main = "Simpson", xlab = "Invertebrate Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyPlant[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[4,], main = "Simpson", xlab = "Plant Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyPlant[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Pielou 
plot(x = PielouInverts, y = RarefiedAlphaDiversityTableOnlyPlant[5,], main = "Pielou's Evenness", xlab = "Invertebrate Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = PielouInverts, y = RarefiedAlphaDiversityTableOnlyPlant[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouPlants, y = RarefiedAlphaDiversityTableOnlyPlant[5,], main = "Pielou's Evenness", xlab = "Plant Diversity", ylab = "Plant Associated Fungi Diversity", pch = 19)
cor.test(x = PielouPlants, y = RarefiedAlphaDiversityTableOnlyPlant[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Shannon 
plot(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[1,], main = "Shannon", xlab = "Invertebrate Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = ShannonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[1,], main = "Shannon", xlab = "Plant Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = ShannonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[1,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Richness 
plot(x = RichnessInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[2,], main = "Richness", xlab = "Invertebrate Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = RichnessInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = RichnessPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[2,], main = "Richness", xlab = "Plant Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = RichnessPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[2,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Chao1 
plot(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[3,], main = "Chao1", xlab = "Invertebrate Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = Chao1DiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[3,], main = "Chao1", xlab = "Plant Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = Chao1DiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[3,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Simpson 
plot(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[4,], main = "Simpson", xlab = "Invertebrate Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = SimpsonDiversityInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[4,], main = "Simpson", xlab = "Plant Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = SimpsonDiversityPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[4,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

#Pielou 
plot(x = PielouInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[5,], main = "Pielou's Evenness", xlab = "Invertebrate Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = PielouInverts, y = RarefiedAlphaDiversityTableOnlyAniPath[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)

plot(x = PielouPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[5,], main = "Pielou's Evenness", xlab = "Plant Diversity", ylab = "Animal Pathogens Diversity", pch = 19)
cor.test(x = PielouPlants, y = RarefiedAlphaDiversityTableOnlyAniPath[5,], alternative = c("two.sided"), method = "spearman", exact = FALSE)


library(ecodist)
#Let's make a subsample for mantel tests
#Now I want to create a single subsample for beta diversity analyses
set.seed(1210)

Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
dim(Subsample)
rowSums(Subsample)
#I should create an Animal Pathogen subsample
APSubsample <- Subsample[,which(FungusKeyReduced[,1]==1)]
#Now a plant associated subsample
PASubsample <- Subsample[,which(FungusKeyReduced[,2]==1 & FungusKeyReduced[,1]==0)]



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

library(VennDiagram)
plot.new()
draw.pairwise.venn(area1 = sum(FungusKey[,1]),
                   area2 = sum(FungusKey[,2]),
                   cross.area = sum(FungusKey[,1] & FungusKey[,2]),
                   category = c("Animal Pathogen", "Plant Associated"),
                   cex = 3, cat.cex = 3, margin =0.30, cat.pos = c(315,135), cat.dist = c(.1,.1), fill = c("red", "green"))

par(mfrow=c(1,1),mar=c(9.1, 9.1, 4.1, 5.1))
