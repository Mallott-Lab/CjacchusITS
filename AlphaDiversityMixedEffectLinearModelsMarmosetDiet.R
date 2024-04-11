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

#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Loading in relevant packages
library(qiime2R)#We need this to bring in the qza objects
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(stringr)#For strings
library(lme4)#For mixed effect models
library(car)#for type III ANOVA
library(agricolae)#for the Tukey HSD

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

#So, in the rarefy vs. subsample debate, these are the definitions that I am using
#Rarefy means to take a subsample several times, calculate alpha diversity, then average those values
#Subsample means to take a set number of reads and treat that as the new sample
#Several of these functions subsample when they should rarefy

#I can write a loop to actually rarefy my data
#I should calculate multiple metrics while I am at it
#This will do Simpson, Shannon, Chao1, Richness, and Pielou's evenness (Shannon divided by ln(Richness))

ShannonTable<- matrix(nrow=1000,ncol=60)
RichnessTable<- matrix(nrow=1000,ncol=60)
ChaoTable<- matrix(nrow=1000,ncol=60)
SimpsonTable<- matrix(nrow=1000,ncol=60)
PielouTable<- matrix(nrow=1000,ncol=60)

set.seed(50)

for (i in 1:nrow(ShannonTable)){
  Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
  ShannonTable[i,]<-diversity(Subsample, index="shannon")
  SimpsonTable[i,]<-diversity(Subsample, index="simpson")
  ChaoTable[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTable[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTable[i,]<- ShannonTable[i,]/log(RichnessTable[i,])
}

RarefiedAlphaDiversityTable <- matrix(nrow = 5, ncol = 60)
rownames(RarefiedAlphaDiversityTable) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTable) <- colnames(FungTableReduced)
RarefiedAlphaDiversityTable[1,] <- colMeans(ShannonTable)
RarefiedAlphaDiversityTable[2,] <- colMeans(RichnessTable)
RarefiedAlphaDiversityTable[3,] <- colMeans(ChaoTable)
RarefiedAlphaDiversityTable[4,] <- colMeans(SimpsonTable)
RarefiedAlphaDiversityTable[5,] <- colMeans(PielouTable)
RarefiedAlphaDiversityTable


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


#First I want to make a data frame with everything in it
ModelData <- data.frame(cbind(MetadataReduced,ShannonDiversityInverts, SimpsonDiversityInverts,Chao1DiversityInverts,RichnessInverts,PielouInverts))
#I want to make sure the factors are coded as factors
ModelData$Sex <-as.factor(ModelData$Sex)
ModelData$Group <-as.factor(ModelData$Group)
ModelData$Season <-as.factor(ModelData$Season)
ModelData$Reproductive <-as.factor(ModelData$Reproductive)
ModelData$Individual <-as.factor(ModelData$Individual)
ModelData$Preservative <-as.factor(ModelData$Preservative)
ModelData$Age <-as.factor(ModelData$Age)
ModelData$Time2 <-as.factor(ModelData$Time2)
#Now I want to put in NAs when appropriate
ModelData$Reproductive[which(ModelData$Reproductive=="Unknown")]<-NA
ModelData$Age[which(ModelData$Age=="Unknown")]<-NA
ModelData$Sex[which(ModelData$Sex=="Unknown")]<-NA

#Now for the mixed effects models
#I am treating season and sex as fixed effects because they are the variables of interest
#Preservative and age are both fixed effects because we sampled all levels, they have few levels
#Group is a random effect because I am not necessarily interested in it directly, but it could matter
#The Time2 (season and year) is nested in season
#Reproductive state could be nested in sex (only females can be pregnant), but this makes it hard to test for significance (can't drop sex while reproductive state is nested in it)
#After dropping the unknown sex individuals, we only have 52 samples, 30 are male (so not-pregnant male) and 10 of the 22 females are "unknown" reproductive state
#I cannot imagine that the information will be helpful, there are a lot of NAs and most reproductive states are poorly sampled, so I am dropping reproductive state completely

#The pattern will be to create a full model, then create a null model that is missing the parameter of interest (one missing sex, the other missing season)
#The full model will be compared to the relevant null model to test for significance of the missing parameter
ShannonModelFull <- lmer(ShannonDiversityInverts ~ 
                           Season +
                           Sex +
                           (1|Group) +
                           #Reproductive +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(ShannonModelFull)
Anova(ShannonModelFull)

ShannonModelNullSeason <- lmer(ShannonDiversityInverts ~ 
                                 #Season +
                                 Sex +
                                 (1|Group) +
                                 #Reproductive +
                                 Preservative +
                                 Age +
                                 Season:Time2,
                               data = ModelData, REML=FALSE)
summary(ShannonModelNullSeason)
ShannonModelNullSex <- lmer(ShannonDiversityInverts ~ 
                              Season +
                              #Sex +
                              (1|Group) +
                              #Reproductive +
                              Preservative +
                              Age +
                              Season:Time2,
                            data = ModelData, REML=FALSE)
summary(ShannonModelNullSex)
#Now for significance tests
anova(ShannonModelNullSeason, ShannonModelFull, test="Chisq")
anova(ShannonModelNullSex, ShannonModelFull, test="Chisq")
#Sex is not significant, Season is not significant


#Same procedure for Richness
RichnessModelFull <- lmer(RichnessInverts ~ 
                            Season +
                            Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(RichnessModelFull)
Anova(RichnessModelFull)

RichnessModelNullSeason <- lmer(RichnessInverts ~ 
                                  #Season +
                                  Sex +
                                  (1|Group) +
                                  #Reproductive +
                                  Preservative +
                                  Age +
                                  Season:Time2,
                                data = ModelData, REML=FALSE)
summary(RichnessModelNullSeason)
RichnessModelNullSex <- lmer(RichnessInverts ~ 
                               Season +
                               #Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(RichnessModelNullSex)
#Tests for significance
anova(RichnessModelNullSeason, RichnessModelFull, test="Chisq")
anova(RichnessModelNullSex, RichnessModelFull, test="Chisq")
#Sex is not significant, Season is not significant

Chao1ModelFull <- lmer(Chao1DiversityInverts ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         #Reproductive +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelData, REML=FALSE)
summary(Chao1ModelFull)
Anova(Chao1ModelFull)

Chao1ModelNullSeason <- lmer(Chao1DiversityInverts ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(Chao1ModelNullSeason)
Chao1ModelNullSex <- lmer(Chao1DiversityInverts ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(Chao1ModelNullSex)
#Test for significance
anova(Chao1ModelNullSeason, Chao1ModelFull, test="Chisq")
anova(Chao1ModelNullSex, Chao1ModelFull, test="Chisq")
#Sex is not significant, season is not significant


SimpsonModelFull <- lmer(SimpsonDiversityInverts ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         #Reproductive +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelData, REML=FALSE)
summary(SimpsonModelFull)
Anova(SimpsonModelFull)

SimpsonModelNullSeason <- lmer(SimpsonDiversityInverts ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(SimpsonModelNullSeason)
SimpsonModelNullSex <- lmer(SimpsonDiversityInverts ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(SimpsonModelNullSex)
#Test for significance
anova(SimpsonModelNullSeason, SimpsonModelFull, test="Chisq")
anova(SimpsonModelNullSex, SimpsonModelFull, test="Chisq")
#Sex is not significant, season is not significant



PielouModelFull <- lm(PielouInverts ~ 
                           Season +
                           Sex +
                           Group +
                           #Reproductive +
                           Preservative +
                           Age,
                         data = ModelData)
summary(PielouModelFull)
#None of the parameters are significant
#The model itself is not significant

#Nothing is significant for invertebrate alpha richness metrics


#Let's look at the plants

ShannonModelFull <- lm(ShannonDiversityPlants ~ 
                           Season +
                           Sex +
                           Group +
                           #Reproductive +
                           Preservative +
                           Age,
                         data = ModelData)
summary(ShannonModelFull)
#None of the parameters are significant
#The overall model is not significant


#Same procedure for Richness
RichnessModelFull <- lmer(RichnessPlants ~ 
                            Season +
                            Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(RichnessModelFull)
Anova(RichnessModelFull)

RichnessModelNullSeason <- lmer(RichnessPlants ~ 
                                  #Season +
                                  Sex +
                                  (1|Group) +
                                  #Reproductive +
                                  Preservative +
                                  Age +
                                  Season:Time2,
                                data = ModelData, REML=FALSE)
summary(RichnessModelNullSeason)
RichnessModelNullSex <- lmer(RichnessPlants ~ 
                               Season +
                               #Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(RichnessModelNullSex)
#Tests for significance
anova(RichnessModelNullSeason, RichnessModelFull, test="Chisq")
anova(RichnessModelNullSex, RichnessModelFull, test="Chisq")
#Sex is ALMOST significant, season is not significant


Chao1ModelFull <- lmer(Chao1DiversityPlants ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         #Reproductive +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelData, REML=FALSE)
summary(Chao1ModelFull)
Anova(Chao1ModelFull)

Chao1ModelNullSeason <- lmer(Chao1DiversityPlants ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(Chao1ModelNullSeason)
Chao1ModelNullSex <- lmer(Chao1DiversityPlants ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(Chao1ModelNullSex)
#Test for significance
anova(Chao1ModelNullSeason, Chao1ModelFull, test="Chisq")
anova(Chao1ModelNullSex, Chao1ModelFull, test="Chisq")
#Sex is not significant, season is not significant


SimpsonModelFull <- lm(SimpsonDiversityPlants ~ 
                           Season +
                           Sex +
                           Group +
                           #Reproductive +
                           Preservative +
                           Age,
                         data = ModelData)
summary(SimpsonModelFull)
#The model is not significant
#None of the parameters are significant

PielouModelFull <- lmer(PielouPlants ~ 
                           Season +
                           Sex +
                           (1|Group) +
                           #Reproductive +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(PielouModelFull)
Anova(PielouModelFull)

PielouModelNullSex <- lmer(PielouPlants ~ 
                          Season +
                          #Sex +
                          (1|Group) +
                          #Reproductive +
                          Preservative +
                          Age +
                          Season:Time2,
                        data = ModelData, REML=FALSE)
summary(PielouModelNullSex)

PielouModelNullSeason <- lmer(PielouPlants ~ 
                          #Season +
                          Sex +
                          (1|Group) +
                          #Reproductive +
                          Preservative +
                          Age +
                          Season:Time2,
                        data = ModelData, REML=FALSE)
summary(PielouModelNullSeason)
#Tests for significance
anova(PielouModelNullSeason, PielouModelFull, test="Chisq")
#Season is not significant
#The sex null model cannot be lmer


PielouModelFull <- lm(PielouPlants ~ 
                          Season +
                          Sex +
                          Group +
                          #Reproductive +
                          Preservative +
                          Age,
                        data = ModelData)
summary(PielouModelFull)
#Preservative is significant, but the model is not