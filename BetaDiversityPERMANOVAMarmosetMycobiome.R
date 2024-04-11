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

#Two of the samples are really throwing the NMDS off
#I am going to save the original subsample then try again without those two samples
#They are rows 43 and 47, one is a pregnant female
SavedSubsample <- Subsample
SavedMetadata <- MetadataReduced

Subsample <- Subsample[-c(43,47),]
dim(Subsample)
dim(MetadataReduced)
MetadataReduced <- MetadataReduced[-c(43,47),]
dim(MetadataReduced)
#I am going to create a Bray-Curtis dissimilarity matrix
BrayDist<-vegdist(Subsample,method="bray",binary=FALSE)
#Here is a Jaccard distance matrix
JaccDist<- vegdist(Subsample,method="jaccard",binary=TRUE)

#Now the data are ready for analyses
#Let's start with the PERMANOVA of the Bray-Curtis dissimilarity matrix
#Dropped Time from the PERMANOVA because it is nested in season and PERMANOVA doesn't quite handle that well
BrayPERMANOVA<-adonis2(BrayDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

BrayPERMANOVA
#Season, preservative, and group are significant in the PERMANOVA

#Let's run NMDS
BrayNMDS<-metaMDS(BrayDist, k = 2,autotransform=FALSE,try=100,trymax=100)

#Let's plot it
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(BrayNMDS$points,pch=19, cex= 0.5,col="white", main= "Bray-Curtis NMDS")
#groups were significant, so let's make those with colors
points(BrayNMDS$points[which(MetadataReduced$Group=="Algaroba"),], pch = 0, cex = 1.5, col = "#E69F00")
points(BrayNMDS$points[which(MetadataReduced$Group=="Coqueiro"),], pch = 1, cex = 1.5, col = "#56B4E9")
points(BrayNMDS$points[which(MetadataReduced$Group=="Cow"),], pch = 2, cex = 1.5, col = "#009E73")
points(BrayNMDS$points[which(MetadataReduced$Group=="F group"),], pch = 15, cex = 1.5, col = "#F0E442")
points(BrayNMDS$points[which(MetadataReduced$Group=="House"),], pch = 17, cex = 1.5, col = "#0072B2")
points(BrayNMDS$points[which(MetadataReduced$Group=="Key"),], pch = 19, cex = 1.5, col = "#D55E00")
points(BrayNMDS$points[which(MetadataReduced$Group=="Princess"),], pch = 8, cex = 1.5, col = "#CC79A7")
points(BrayNMDS$points[which(MetadataReduced$Group=="Road"),], pch = 4, cex = 1.5, col = "#999999")
ordiellipse(BrayNMDS$points,MetadataReduced$Group,col=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"))

plot.new()
legend("center",c("Algaroba","Coqueiro", "Cow", "F group", "House","Key","Princess","Road"),
       pch=c(0,1,2,15,17,19,8,4)
       ,col=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
       ,xjust=0,yjust=0,title = NA,ncol = 2,pt.cex = 1.5)

#Let's plot Preservative
plot(BrayNMDS$points,pch=19, cex= 0.5, main= "Bray-Curtis NMDS")
#groups were significant, so let's mark those with colors
points(BrayNMDS$points[which(MetadataReduced$Preservative=="Ethanol"),], pch = 19, cex = 1.5, col = "#E69F00")
points(BrayNMDS$points[which(MetadataReduced$Preservative=="RNAlater"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(BrayNMDS$points,MetadataReduced$Preservative,col=c("#E69F00","#56B4E9"))


legend("topright",c("Ethanol","RNAlater")
       ,fill=c("#E69F00","#56B4E9")
       ,xjust=0,yjust=0,title = "Preservative",ncol = 1,pt.cex = 1.5)

#Let's plot Season, it approaches significance
plot(BrayNMDS$points,pch=19, cex= 0.5, main= "Bray-Curtis NMDS")
#groups were significant, so let's mark those with colors
points(BrayNMDS$points[which(MetadataReduced$Season=="Dry"),], pch = 19, cex = 1.5, col = "#E69F00")
points(BrayNMDS$points[which(MetadataReduced$Season=="Wet"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(BrayNMDS$points,MetadataReduced$Season,col=c("#E69F00","#56B4E9"))

legend("topright",c("Dry","Wet")
       ,fill=c("#E69F00","#56B4E9")
       ,xjust=0,yjust=0,title = "Season",ncol = 1,pt.cex = 1.5)

#I want to look at betadisper for season
SeasonBrayBetaDisper <- betadisper(d=BrayDist,group=MetadataReduced$Season)
anova(SeasonBrayBetaDisper)
TukeyHSD(SeasonBrayBetaDisper)
#The dispersion of points is not significantly different

#Let's move on to the Jaccard distance PERMANOVA
JaccPERMANOVA<-adonis2(JaccDist~ MetadataReduced$Age + MetadataReduced$Sex + MetadataReduced$Season + MetadataReduced$Group + MetadataReduced$Preservative,by="margin")

JaccPERMANOVA
#Season, Group, and Preservative are significant

#Let's run NMDS
JaccNMDS<-metaMDS(JaccDist, k = 2,autotransform=FALSE,try=100,trymax=100)

#Let's plot it
plot(JaccNMDS$points,pch=19, cex= 0.5,col="white", main= "Jaccard NMDS")
#groups were significant, so let's makr those with colors
points(JaccNMDS$points[which(MetadataReduced$Group=="Algaroba"),], pch = 0, cex = 1.5, col = "#E69F00")
points(JaccNMDS$points[which(MetadataReduced$Group=="Coqueiro"),], pch = 1, cex = 1.5, col = "#56B4E9")
points(JaccNMDS$points[which(MetadataReduced$Group=="Cow"),], pch = 2, cex = 1.5, col = "#009E73")
points(JaccNMDS$points[which(MetadataReduced$Group=="F group"),], pch = 15, cex = 1.5, col = "#F0E442")
points(JaccNMDS$points[which(MetadataReduced$Group=="House"),], pch = 17, cex = 1.5, col = "#0072B2")
points(JaccNMDS$points[which(MetadataReduced$Group=="Key"),], pch = 19, cex = 1.5, col = "#D55E00")
points(JaccNMDS$points[which(MetadataReduced$Group=="Princess"),], pch = 8, cex = 1.5, col = "#CC79A7")
points(JaccNMDS$points[which(MetadataReduced$Group=="Road"),], pch = 4, cex = 1.5, col = "#999999")
ordiellipse(JaccNMDS$points,MetadataReduced$Group,col=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"))

plot.new()
legend("center",c("Algaroba","Coqueiro", "Cow", "F group", "House","Key","Princess","Road"),
       pch=c(0,1,2,15,17,19,8,4)
       ,col=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
       ,xjust=0,yjust=0,title = NA,ncol = 2,pt.cex = 1.5)


#Let's plot Season
plot(JaccNMDS$points,pch=19, cex= 0.5, main= "Jaccard NMDS")
#groups were significant, so let's mark those with colors
points(JaccNMDS$points[which(MetadataReduced$Season=="Dry"),], pch = 19, cex = 1.5, col = "#E69F00")
points(JaccNMDS$points[which(MetadataReduced$Season=="Wet"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(JaccNMDS$points,MetadataReduced$Season,col=c("#E69F00","#56B4E9"))

legend("bottomright",c("Dry","Wet")
       ,fill=c("#E69F00","#56B4E9")
       ,xjust=0,yjust=0,title = "Season",ncol = 1,pt.cex = 1.5)

SeasonJaccBetaDisper <- betadisper(d=JaccDist,group=MetadataReduced$Season)
anova(SeasonJaccBetaDisper)
TukeyHSD(SeasonJaccBetaDisper)
#The dispersion of points is not significantly different

#Let's plot Preservative
plot(JaccNMDS$points,pch=19, cex= 0.5, main= "Jaccard NMDS")
#groups were significant, so let's mark those with colors
points(JaccNMDS$points[which(MetadataReduced$Preservative=="Ethanol"),], pch = 19, cex = 1.5, col = "#E69F00")
points(JaccNMDS$points[which(MetadataReduced$Preservative=="RNAlater"),], pch = 19, cex = 1.5, col = "#56B4E9")
ordiellipse(JaccNMDS$points,MetadataReduced$Preservative,col=c("#E69F00","#56B4E9"))

legend("bottomright",c("Ethanol","RNAlater")
       ,fill=c("#E69F00","#56B4E9")
       ,xjust=0,yjust=0,title = "Preservative",ncol = 1,pt.cex = 1.5)

