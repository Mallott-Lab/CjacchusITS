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

#I should figure out which ASVs remain in the analyses









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

#I want to drop things that are no longer in the FungTableReduced
FungusKeyReduced <- FungusKey[-which(rowSums(FungTableReduced)==0),]

FungTableReduced <- FungTableReduced[-which(rowSums(FungTableReduced)==0),]
dim(FungTableReduced)
#We lose 34 plant pathogens and 59 plant associated when we drop samples etc.
length(which(rowSums(FungusKeyReduced)==2))
#There are 189 remaining ASVs that are both an animal pathogen and plant associated
colSums(FungusKeyReduced)
#We have 244 Animal pathogens and 397 plant associated


#I want to calculate Alpha diversity metrics associated with Animal Pathogens and Plant Associated Fungi
#There are ASVs that are both animal pathogens AND plant associated
#I am going to run the analyses in three different ways
#First way, I will include ASVs that are both as Animal Pathogens
#Second way, I will include ASVs that are both as Plant Associated
#Third way, I will exclude ASVs that are both


ShannonTable<- matrix(nrow=1000,ncol=52)
RichnessTable<- matrix(nrow=1000,ncol=52)
ChaoTable<- matrix(nrow=1000,ncol=52)
SimpsonTable<- matrix(nrow=1000,ncol=52)
PielouTable<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTable)){
  Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,1]==1)]
  ShannonTable[i,]<-diversity(Subsample, index="shannon")
  SimpsonTable[i,]<-diversity(Subsample, index="simpson")
  ChaoTable[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTable[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTable[i,]<- ShannonTable[i,]/log(RichnessTable[i,])
}


RarefiedAlphaDiversityTableBothAsAniPath <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableBothAsAniPath) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableBothAsAniPath) <- colnames(FungTableReduced)
RarefiedAlphaDiversityTableBothAsAniPath[1,] <- colMeans(ShannonTable)
RarefiedAlphaDiversityTableBothAsAniPath[2,] <- colMeans(RichnessTable)
RarefiedAlphaDiversityTableBothAsAniPath[3,] <- colMeans(ChaoTable)
RarefiedAlphaDiversityTableBothAsAniPath[4,] <- colMeans(SimpsonTable)
RarefiedAlphaDiversityTableBothAsAniPath[5,] <- colMeans(PielouTable)
RarefiedAlphaDiversityTableBothAsAniPath


#Now I want to count them both as Plant Associated

ShannonTable<- matrix(nrow=1000,ncol=52)
RichnessTable<- matrix(nrow=1000,ncol=52)
ChaoTable<- matrix(nrow=1000,ncol=52)
SimpsonTable<- matrix(nrow=1000,ncol=52)
PielouTable<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTable)){
  Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,2]==1)]
  ShannonTable[i,]<-diversity(Subsample, index="shannon")
  SimpsonTable[i,]<-diversity(Subsample, index="simpson")
  ChaoTable[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTable[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTable[i,]<- ShannonTable[i,]/log(RichnessTable[i,])
}


RarefiedAlphaDiversityTableBothAsPlant <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableBothAsPlant) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableBothAsPlant) <- colnames(FungTableReduced)
RarefiedAlphaDiversityTableBothAsPlant[1,] <- colMeans(ShannonTable)
RarefiedAlphaDiversityTableBothAsPlant[2,] <- colMeans(RichnessTable)
RarefiedAlphaDiversityTableBothAsPlant[3,] <- colMeans(ChaoTable)
RarefiedAlphaDiversityTableBothAsPlant[4,] <- colMeans(SimpsonTable)
RarefiedAlphaDiversityTableBothAsPlant[5,] <- colMeans(PielouTable)
RarefiedAlphaDiversityTableBothAsPlant


#Now I want to only count the ones that are just Plant Associated and not also Animal Pathogens

ShannonTable<- matrix(nrow=1000,ncol=52)
RichnessTable<- matrix(nrow=1000,ncol=52)
ChaoTable<- matrix(nrow=1000,ncol=52)
SimpsonTable<- matrix(nrow=1000,ncol=52)
PielouTable<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTable)){
  Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,2]==1 & FungusKeyReduced[,1]==0)]
  ShannonTable[i,]<-diversity(Subsample, index="shannon")
  SimpsonTable[i,]<-diversity(Subsample, index="simpson")
  ChaoTable[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTable[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTable[i,]<- ShannonTable[i,]/log(RichnessTable[i,])
}


RarefiedAlphaDiversityTableOnlyPlant <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableOnlyPlant) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableOnlyPlant) <- colnames(FungTableReduced)
RarefiedAlphaDiversityTableOnlyPlant[1,] <- colMeans(ShannonTable)
RarefiedAlphaDiversityTableOnlyPlant[2,] <- colMeans(RichnessTable)
RarefiedAlphaDiversityTableOnlyPlant[3,] <- colMeans(ChaoTable)
RarefiedAlphaDiversityTableOnlyPlant[4,] <- colMeans(SimpsonTable)
RarefiedAlphaDiversityTableOnlyPlant[5,] <- colMeans(PielouTable)
RarefiedAlphaDiversityTableOnlyPlant


#Now I want to only count the ones that are just Animal Pathogens and not also Plant Associated


ShannonTable<- matrix(nrow=1000,ncol=52)
RichnessTable<- matrix(nrow=1000,ncol=52)
ChaoTable<- matrix(nrow=1000,ncol=52)
SimpsonTable<- matrix(nrow=1000,ncol=52)
PielouTable<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTable)){
  Subsample<-rrarefy(t(FungTableReduced),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,2]==0 & FungusKeyReduced[,1]==1)]
  ShannonTable[i,]<-diversity(Subsample, index="shannon")
  SimpsonTable[i,]<-diversity(Subsample, index="simpson")
  ChaoTable[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTable[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTable[i,]<- ShannonTable[i,]/log(RichnessTable[i,])
}


RarefiedAlphaDiversityTableOnlyAniPath <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableOnlyAniPath) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableOnlyAniPath) <- colnames(FungTableReduced)
RarefiedAlphaDiversityTableOnlyAniPath[1,] <- colMeans(ShannonTable)
RarefiedAlphaDiversityTableOnlyAniPath[2,] <- colMeans(RichnessTable)
RarefiedAlphaDiversityTableOnlyAniPath[3,] <- colMeans(ChaoTable)
RarefiedAlphaDiversityTableOnlyAniPath[4,] <- colMeans(SimpsonTable)
RarefiedAlphaDiversityTableOnlyAniPath[5,] <- colMeans(PielouTable)
RarefiedAlphaDiversityTableOnlyAniPath




#Now I need to see what predicts Animal Pathogen Alpha Diversity when I don't count the ASVs that are both animal and plant associated
ModelData <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableOnlyAniPath)))
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




ShannonModelFull <- lm(Shannon ~ 
                           Season +
                           Sex +
                           Group +
                           #Reproductive +
                           Preservative +
                           Age,
                         data = ModelData)
summary(ShannonModelFull)
#The model itself is not significant (p-value = 0.4309)
#Sex is the only significant parameter, with male being different from female

T3ANOVAShannon <- Anova(ShannonModelFull , type = 3)

TukeyHSDT3ANOVAShannon <- HSD.test(ShannonModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)

View(TukeyHSDT3ANOVAShannon$comparison)
#None of the comparisons are significant


RichnessModelFull <- lm(Richness ~ 
                            Season +
                            Sex +
                            Group +
                            #Reproductive +
                            Preservative +
                            Age,
                          data = ModelData)
summary(RichnessModelFull)
#The model itself is not significant (p-value = 0.882)
#No parameter is significant



Chao1ModelFull <- lm(Chao1 ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(Chao1ModelFull)
#The model itself is not significant (p-value = 0.886)
#No parameter is significant

SimpsonModelFull <- lm(Simpson ~ 
                            Season +
                            Sex +
                            Group +
                            #Reproductive +
                            Preservative +
                            Age,
                          data = ModelData)
summary(SimpsonModelFull)
#The model itself is not significant (p-value = 0.6564)
#No parameter is significant

PielouModelFull <- lm(Pielou ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(PielouModelFull)
#The model itself is not significant (p-value = 0.5152)
#The Road and F group groups are different from Abocarpa or whatever it is called

T3ANOVAPielou <- Anova(PielouModelFull , type = 3)

TukeyHSDT3ANOVAPielou <- HSD.test(PielouModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)

View(TukeyHSDT3ANOVAPielou$comparison)




#Let's see if it changes when I look at all Animal Pathogens
ModelData <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableBothAsAniPath)))
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

#Shannon
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
RichnessModelFull <- lm(Richness ~ 
                            Season +
                            Sex +
                            Group +
                            #Reproductive +
                            Preservative +
                            Age,
                          data = ModelData)
summary(RichnessModelFull)

#Running a type 3 ANOVA to get a better idea
T3ANOVARichness <- Anova(RichnessModelFull , type = 3)

TukeyHSDT3ANOVARichness <- HSD.test(RichnessModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)
View(TukeyHSDT3ANOVARichness$comparison)


Chao1ModelFull <- lm(Chao1 ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(Chao1ModelFull)
#Running a type 3 ANOVA to get a better idea
T3ANOVAChao1 <- Anova(Chao1ModelFull , type = 3)

TukeyHSDT3ANOVAChao <- HSD.test(Chao1ModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)
View(TukeyHSDT3ANOVAChao$comparison)


SimpsonModelFull <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(SimpsonModelFull)
T3ANOVASimpson <- Anova(SimpsonModelFull , type = 3)

TukeyHSDT3ANOVASimpson <- HSD.test(SimpsonModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)
View(TukeyHSDT3ANOVASimpson$comparison)


PielouModelFull <- lm(Pielou ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(PielouModelFull)


#Next I need to see how the various ways of considering the plants work for them


#Now I need to see what predicts Plant Associated Alpha Diversity when I don't count the ASVs that are both animal and plant associated
ModelData <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableOnlyPlant)))
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


#Shannon
ShannonModelFull <- lm(Shannon ~ 
                           Season +
                           Sex +
                           Group +
                           #Reproductive +
                           Preservative +
                           Age,
                         data = ModelData)
summary(ShannonModelFull)

#Richness
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
#Now for significance tests
anova(RichnessModelNullSeason, RichnessModelFull, test="Chisq")
anova(RichnessModelNullSex, RichnessModelFull, test="Chisq")
#Sex is significant, season is not




#Chao1
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
#Now for significance tests
anova(Chao1ModelNullSeason, Chao1ModelFull, test="Chisq")
anova(Chao1ModelNullSex, Chao1ModelFull, test="Chisq")
#Nothing is significant



#Simpson
SimpsonModelFull <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(SimpsonModelFull)

#Pielou
PielouModelFull <- lm(Pielou ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(PielouModelFull)
T3ANOVAPielou <- Anova(PielouModelFull , type = 3)

TukeyHSDT3ANOVAPielou <- HSD.test(PielouModelFull, c("Season","Sex", "Group", "Preservative", "Age"), group = FALSE)
View(TukeyHSDT3ANOVAPielou$comparison)

#Now I need to see what predicts Plant Associated Alpha Diversity when I count the ASVs that are both animal and plant associated as plant associated
ModelData <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableBothAsPlant)))
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

#Shannon
ShannonModelFull <- lm(Shannon ~ 
                            Season +
                            Sex +
                            Group +
                            #Reproductive +
                            Preservative +
                            Age,
                          data = ModelData)
summary(ShannonModelFull)

#Richness
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
#Now for significance tests
anova(RichnessModelNullSeason, RichnessModelFull, test="Chisq")
anova(RichnessModelNullSex, RichnessModelFull, test="Chisq")
#Sex is significant, season is not

#Chao1
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
#Now for significance tests
anova(Chao1ModelNullSeason, Chao1ModelFull, test="Chisq")
anova(Chao1ModelNullSex, Chao1ModelFull, test="Chisq")
#Sex is significant, season is not


#Simpson
SimpsonModelFull <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(SimpsonModelFull)

#Pielou
PielouModelFull <- lm(Pielou ~ 
                         Season +
                         Sex +
                         Group +
                         #Reproductive +
                         Preservative +
                         Age,
                       data = ModelData)
summary(PielouModelFull)
