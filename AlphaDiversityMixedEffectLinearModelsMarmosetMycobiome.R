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

#Now it is time to look at the effects of my parameters on various alpha diversity metrics

#First I want to make a data frame with everything in it
RarefiedAlphaDiversityTable
MetadataReduced
ModelData <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTable)))
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
#I need to drop samples with unknown sex, it is a major question here
ModelData<-ModelData[-which(is.na(ModelData$Sex)==TRUE),]


#Now for the mixed effects models
#I am treating season and sex as fixed effects because they are the variables of interest
#Preservative and age are both fixed effects because we sampled all levels, they have few levels
#Group is a random effect because I am not necessarily interested in it directly, but it could matter
#The Time2 (season and year) is nested in season
#Reproductive state could be nested in sex (only females can be pregnant), but this makes it hard to test for significance (can't drop sex while reproductive state is nested in it)
#After dropping the unknown sex individuals, we only have 52 samples, 30 are male (so not-pregnant male) and 10 of the 22 females are "unknown" reproductive state
#I cannot imagine that the information will be helpful, there are a lot of NAs and most reproductive states are poorly stampled, so I am dropping reproductive state completely

#The pattern will be to create a full model, then create a null model that is missing the parameter of interest (one missing sex, the other missing season)
#The full model will be compared to the relevant null model to test for significance of the missing parameter
ShannonModelFull <- lmer(Shannon ~ 
                       Season +
                       Sex +
                       (1|Group) +
                       #Reproductive +
                       Preservative +
                       Age +
                       Season:Time2,
                     data = ModelData, REML=FALSE)
summary(ShannonModelFull)


ShannonModelNullSeason <- lmer(Shannon ~ 
                           #Season +
                           Sex +
                           (1|Group) +
                           #Reproductive +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(ShannonModelNullSeason)
ShannonModelNullSex <- lmer(Shannon ~ 
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
#Sex is significant, season is not

#Same procedure for Richness
RichnessModelFull <- lmer(Richness ~ 
                           Season +
                           Sex +
                           (1|Group) +
                           #Reproductive +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(RichnessModelFull)


RichnessModelNullSeason <- lmer(Richness ~ 
                                 #Season +
                                 Sex +
                                 (1|Group) +
                                 #Reproductive +
                                 Preservative +
                                 Age +
                                 Season:Time2,
                               data = ModelData, REML=FALSE)
summary(RichnessModelNullSeason)
RichnessModelNullSex <- lmer(Richness ~ 
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
#Sex is significant, season is not

Chao1ModelFull <- lmer(Chao1 ~ 
                            Season +
                            Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(Chao1ModelFull)


Chao1ModelNullSeason <- lmer(Chao1 ~ 
                                  #Season +
                                  Sex +
                                  (1|Group) +
                                  #Reproductive +
                                  Preservative +
                                  Age +
                                  Season:Time2,
                                data = ModelData, REML=FALSE)
summary(Chao1ModelNullSeason)
Chao1ModelNullSex <- lmer(Chao1 ~ 
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
#Sex is significant, season is not


SimpsonModelFull <- lmer(Simpson ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         #Reproductive +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelData, REML=FALSE)
summary(SimpsonModelFull)


SimpsonModelNullSeason <- lmer(Simpson ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               #Reproductive +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(SimpsonModelNullSeason)
SimpsonModelNullSex <- lmer(Simpson ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            #Reproductive +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(SimpsonModelNullSex)
#Now to test for significance
anova(SimpsonModelNullSeason, SimpsonModelFull, test="Chisq")
anova(SimpsonModelNullSex, SimpsonModelFull, test="Chisq")
#Sex is significant, season is not

PielouModelFull <- lm(Pielou ~ 
                           Season +
                           Sex +
                           #(1|Group) +
                           #Reproductive +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(PielouModelFull)


PielouModelNullSeason <- lmer(Pielou ~ 
                                 #Season +
                                 Sex +
                                 #(1|Group) +
                                 #Reproductive +
                                 Preservative +
                                 Age +
                                 Season:Time2,
                               data = ModelData, REML=FALSE)
summary(PielouModelNullSeason)
PielouModelNullSex <- lmer(Pielou ~ 
                              Season +
                              #Sex +
                              #(1|Group) +
                              #Reproductive +
                              Preservative +
                              Age +
                              Season:Time2,
                            data = ModelData, REML=FALSE)
summary(PielouModelNullSex)
#Now to test for significance
anova(PielouModelNullSeason, PielouModelFull, test="Chisq")
anova(PielouModelNullSex, PielouModelFull, test="Chisq")
#Sex and season are both not significant


#I want to visualize these results, probably with box plots for each alpha diversity metric


par(mfrow=c(2,2),mar=c(5.1, 5.1, 4.1, 2.1))
boxplot(
  ModelData$Shannon[which(ModelData$Sex=="Male")],
  ModelData$Shannon[which(ModelData$Sex=="Female")]
  ,xaxt="n",ylab="Mycobiome Shannon Diversity",xlab="Sex",
  cex.lab =1.8, cex=3,ylim=c(0,5),lwd=2)
axis(side = 1, at=c(1,2),labels=c("Male","Female"))

boxplot(
  ModelData$Richness[which(ModelData$Sex=="Male")],
  ModelData$Richness[which(ModelData$Sex=="Female")]
  ,xaxt="n",ylab="Mycobiome Richness Diversity",xlab="Sex",
  cex.lab =1.8, cex=3,ylim=c(0,150),lwd=2)
axis(side = 1, at=c(1,2),labels=c("Male","Female"))

boxplot(
  ModelData$Chao1[which(ModelData$Sex=="Male")],
  ModelData$Chao1[which(ModelData$Sex=="Female")]
  ,xaxt="n",ylab="Mycobiome Chao1 Diversity",xlab="Sex",
  cex.lab =1.8, cex=3,ylim=c(0,150),lwd=2)
axis(side = 1, at=c(1,2),labels=c("Male","Female"))

boxplot(
  ModelData$Simpson[which(ModelData$Sex=="Male")],
  ModelData$Simpson[which(ModelData$Sex=="Female")]
  ,xaxt="n",ylab="Mycobiome Simpson Diversity",xlab="Sex",
  cex.lab =1.8, cex=3,ylim=c(0,1.1),lwd=2)
axis(side = 1, at=c(1,2),labels=c("Male","Female"))

par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
boxplot(
  ModelData$Pielou[which(ModelData$Sex=="Male")],
  ModelData$Pielou[which(ModelData$Sex=="Female")]
  ,xaxt="n",ylab="Mycobiome Pielou's Evenness",xlab="Sex",
  cex.lab =1.8, cex=3,ylim=c(0,1.1),lwd=2)
axis(side = 1, at=c(1,2),labels=c("Male","Female"))
