# Loading packages----
#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs")
#Loading in relevant packages
library(qiime2R)#We need this to bring in the qza objects
library(vegan)#For various indices and metrics
library(fossil)#For Chao1
library(stringr)#For working with strings
library(lme4)#For mixed effect models
library(car)#for type III ANOVA
library(ggnewscale)#For plots
library(ANCOMBC)#For ANCOM-BC
library(phyloseq)#For ANCOM-BC
library(ggtext)#For Plots
library(ggplot2)#For Plots
library(glue)#For Plots
library(cowplot)#For Plots
library(ggordiplots)#For Plots

# Loading files----
#Let's bring in the objects from Qiime2
#These are mycobiome tables
tabletrimmedQiime2Object<-read_qza("table-trimmed.qza", rm = FALSE)
taxtableQiime2Object<-read_qza("taxonomy.qza")

#Let's pull out what I specifically want to work with
FungTable<-tabletrimmedQiime2Object$data
TaxTable<-taxtableQiime2Object$data
rownames(TaxTable) <- TaxTable$Feature.ID
#I output these files to run things through FUNGuild
#write.table(TaxTable,"TaxTable.txt", sep = "\t")
#I had to delete the column of index numbers in this
#write.table(FungTable,"ASVTable.txt", sep = "\t")

#Now I want to bring in my metadata
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork")
metadata<- read.table("Cjacchus_ITS_metadata.txt",sep="\t",header=TRUE)

# Formatting metadata----
#This metadata file is not in the same sequence as the table
metadata<- metadata[order(metadata$sampleid),]
#I want to make a change to reproductive status to distinguish non-pregnant females from males
metadata$Reproductive[which(metadata$Sex=="Male")]<-"NPMale"
table(metadata$Reproductive)
#I want to drop empty columns
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

# Calculating within sample diversity metrics of the overall mycobiome----
#So, in the rarefy vs. subsample debate, these are the definitions that I am using
#Rarefy means to take a subsample several times, calculate alpha diversity, then average those values
#Subsample means to take a set number of reads and treat that as the new sample
#Several of these functions subsample when they should rarefy

#I can write a loop to actually rarefy my data
#I should calculate multiple metrics while I am at it
#This will do Simpson, Shannon, Chao1, Richness, and Pielou's evenness (Shannon divided by ln(Richness))
#These will inform within sample diversity analyses of the marmoset mycobiome

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

# Fixed and Mixed Effects Linear Regressions Testing Correlations between sex/season and within sample diversity metrics----
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
#I need to drop samples with unknown sex, it is a major parameter of interest here
ModelData<-ModelData[-which(is.na(ModelData$Sex)==TRUE),]

#This code will create a full model, then create a null model that is missing the parameter of interest (one missing sex, the other missing season)
#The full model will be compared to the relevant null model to test for significance of the missing parameter
#This will test the Shannon diversity of the overall mycobiome
ShannonModelFull <- lmer(Shannon ~ 
                           Season +
                           Sex +
                           (1|Group) +
                           Preservative +
                           Age +
                           Season:Time2,
                         data = ModelData, REML=FALSE)
summary(ShannonModelFull)
Anova(ShannonModelFull)

ShannonModelNullSeason <- lmer(Shannon ~ 
                                 #Season +
                                 Sex +
                                 (1|Group) +
                                 Preservative +
                                 Age +
                                 Season:Time2,
                               data = ModelData, REML=FALSE)
summary(ShannonModelNullSeason)
ShannonModelNullSex <- lmer(Shannon ~ 
                              Season +
                              #Sex +
                              (1|Group) +
                              Preservative +
                              Age +
                              Season:Time2,
                            data = ModelData, REML=FALSE)
summary(ShannonModelNullSex)
#Now for significance tests
anova(ShannonModelNullSeason, ShannonModelFull, test="Chisq")
anova(ShannonModelNullSex, ShannonModelFull, test="Chisq")
#Sex is significant, season is not

#Here are box and whiskers plots showing these results
PlotA <- ggplot(ModelData, aes(x = Season, y = Shannon, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Shannon diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  annotate("text",x=0.5,y=3.9,label="italic(p) - value == 0.173", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.5,y=3.6,label= expression( chi ^"2"~ "= 1.858"), hjust =0, size=6) +
  annotate("text",x=0.5,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  theme(aspect.ratio = 1) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotA

PlotA2 <- ggplot(ModelData, aes(x = Sex, y = Shannon, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Shannon diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  annotate("text",x=0.5,y=3.9,label="italic(p) - value < 0.001", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.5,y=3.6,label= expression( chi ^"2"~ "= 12.077"), hjust =0, size=6) +
  annotate("text",x=0.5,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE) +
  geom_segment(aes(x=1,y=4,xend=2,yend=4),inherit.aes = FALSE) +
  theme(aspect.ratio = 1) +
  annotate("text",x=1.5,y=4.1,label="***", size = 10)
PlotA2

#Same procedure for overall mycobiome Richness
RichnessModelFull <- lmer(Richness ~ 
                            Season +
                            Sex +
                            (1|Group) +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(RichnessModelFull)
Anova(RichnessModelFull)

RichnessModelNullSeason <- lmer(Richness ~ 
                                  #Season +
                                  Sex +
                                  (1|Group) +
                                  Preservative +
                                  Age +
                                  Season:Time2,
                                data = ModelData, REML=FALSE)
summary(RichnessModelNullSeason)
RichnessModelNullSex <- lmer(Richness ~ 
                               Season +
                               #Sex +
                               (1|Group) +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(RichnessModelNullSex)
#Tests for significance
anova(RichnessModelNullSeason, RichnessModelFull, test="Chisq")
anova(RichnessModelNullSex, RichnessModelFull, test="Chisq")
#Sex is significant, season is not

PlotB <- ggplot(ModelData, aes(x = Season, y = Richness, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("ASV Richness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  theme(legend.position="bottom") +
  annotate("text",x=0.05,y=150,label="italic(p) - value == 0.292", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=135,label= expression( chi ^"2"~ "= 1.108"), hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  theme(aspect.ratio = 1) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotB

#This creates a legend that will be used later when generating Figure 1
PlotBLegend <- ggplot(ModelData, aes(x = Season, y = Richness, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  ylab("ASV Richness")+
  labs(color = "Climatic Period") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  theme(legend.position="bottom") +
  annotate("text",x=0.05,y=150,label="italic(p) - value == 0.292", parse=TRUE, hjust =0, size=5) +
  annotate("text",x=0.05,y=140,label= expression( chi ^"2"~ "= 1.1075"), hjust =0, size=5) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotBLegend


PlotB2 <- ggplot(ModelData, aes(x = Sex, y = Richness, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("ASV Richness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  theme(legend.position = "bottom") +
  geom_segment(aes(x=1,y=145,xend=2,yend=145),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=150,label="**", size = 10) +
  annotate("text",x=0.05,y=170,label="italic(p) - value == 0.009", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=158,label=expression( chi ^"2" ~ "= 6.828"), hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  theme(aspect.ratio = 1) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotB2

PlotB2Legend <- ggplot(ModelData, aes(x = Sex, y = Richness, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  ylab("ASV Richness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  theme(legend.position = "bottom") +
  geom_segment(aes(x=1,y=145,xend=2,yend=145),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=150,label="**", size = 10) +
  annotate("text",x=0.05,y=160,label="italic(p) - value == 0.009", parse=TRUE, hjust =0, size=5) +
  annotate("text",x=0.05,y=152,label=expression( chi ^"2" ~ "= 6.8275"), hjust =0, size=5) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotB2Legend

#Same procedure for Chao1
Chao1ModelFull <- lmer(Chao1 ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelData, REML=FALSE)
summary(Chao1ModelFull)
Anova(Chao1ModelFull)

Chao1ModelNullSeason <- lmer(Chao1 ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelData, REML=FALSE)
summary(Chao1ModelNullSeason)
Chao1ModelNullSex <- lmer(Chao1 ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelData, REML=FALSE)
summary(Chao1ModelNullSex)
#Test for significance
anova(Chao1ModelNullSeason, Chao1ModelFull, test="Chisq")
anova(Chao1ModelNullSex, Chao1ModelFull, test="Chisq")
#Sex is significant, season is not

PlotSup1 <- ggplot(ModelData, aes(x = Sex, y = Chao1, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Chao1 Diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  geom_segment(aes(x=1,y=150,xend=2,yend=150),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=155,label="**", size = 10) +
  annotate("text",x=0.25,y=160,label="italic(p) - value == 0.009", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.25,y=152,label=expression( chi ^"2" ~ "= 6.8218"), hjust =0, size=6) +
  annotate("text",x=0.2,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotSup1

PlotSup2 <- ggplot(ModelData, aes(x = Season, y = Chao1, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Chao1 Diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  annotate("text",x=0.25,y=160,label="italic(p) - value == 0.233", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.25,y=152,label=expression( chi ^"2" ~ "= 1.424"), hjust =0, size=6) +
  annotate("text",x=0.2,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotSup2


#When I attempted to run the same process as the previous three metrics the model was singular
#The Group, when treated as a random effect, was creating the problem
#I made group into a fixed effect, which changed the model from a mixed effects model to a fixed effects linear model (hence the change in function)
#This model had issues with aliasing, so I had to drop Time2 as a parameter too
SimpsonModelFull <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         Preservative +
                         Age,
                       data = ModelData)
summary(SimpsonModelFull)


#The overall model itself is not significant (p = 0.1159)
#Wet season is significantly different from dry season (p = 0.03851)
#Male is significantly different from female (p = 0.00772)
#Subadult is significantly different from adult (p = 0.03398)

#Plots of the Simpson diversity index
PlotC <- ggplot(ModelData, aes(x = Season, y = Simpson, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Simpson diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  geom_segment(aes(x=1,y=1,xend=2,yend=1),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=1.03,label="*", size = 10) +
  annotate("text",x=0.05,y=1.25,label="italic(p) - value == 0.039", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=1.17,label="italic(t) - value == -2.144", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  theme(aspect.ratio = 1) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotC

PlotC2 <- ggplot(ModelData, aes(x = Sex, y = Simpson, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Simpson diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  geom_segment(aes(x=1,y=1,xend=2,yend=1),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=1.03,label="**", size = 10) +
  annotate("text",x=0.05,y=1.25,label="italic(p) - value == 0.008", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=1.17,label="italic(t) - value == 2.814", parse=TRUE,hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  theme(aspect.ratio = 1) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotC2


#Now for Pielou's evenness, which also uses a fixed effects approach
PielouModelFull <- lm(Pielou ~ 
                        Season +
                        Sex +
                        Group +
                        Preservative +
                        Age,
                      data = ModelData)
summary(PielouModelFull)
#The overall model is not significant (p = 0.1331)
#Male is still significantly different from female (p = 0.00998)
#Subadult is still significantly different from Adult (p = 0.03134)

PlotD <- ggplot(ModelData, aes(x = Season, y = Pielou, color = Season)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Pielou's Evenness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  annotate("text",x=0.05,y=1.1,label="italic(p) - value == 0.099", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=1.02,label="italic(t) - value == -1.693", parse=TRUE,hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotD

PlotD2 <- ggplot(ModelData, aes(x = Sex, y = Pielou, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  ylab("Pielou's Evenness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot(show.legend = FALSE) +
  geom_segment(aes(x=1,y=0.9,xend=2,yend=0.9),inherit.aes = FALSE) +
  annotate("text",x=1.5,y=0.93,label="**", size = 10) +
  annotate("text",x=0.05,y=1.1,label="italic(p) - value == 0.01", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=0.05,y=1.02,label="italic(t) - value == 2.712", parse=TRUE,hjust =0, size=6) +
  annotate("text",x=0,y=1.05,label="", parse=TRUE, hjust =0, size = 10) +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotD2

#I am creating legends for Figure 1
legendSeason <- ggdraw(cowplot::get_plot_component(PlotBLegend, 'guide-box-bottom', return_all = TRUE))
legendSex <- ggdraw(cowplot::get_plot_component(PlotB2Legend, 'guide-box-bottom', return_all = TRUE))

#Supplemental figure
plot_grid(PlotD2 + theme(legend.position="none"),
          PlotD + theme(legend.position="none"),
          legendSex,legendSeason,nrow=2,ncol=2,
          rel_heights = c(1,.3),
          labels = c("A","B"), align = "hv")

plot_grid(PlotSup1 + theme(legend.position="none"),
          PlotSup2 + theme(legend.position="none"),
          PlotA2,PlotA,legendSex,legendSeason,nrow=3,ncol=2,
          rel_heights = c(1,1,.3),
          labels = c("A","B","C","D"), align = "hv")

# Calculating between sample diversity metrics for the overall mycobiome----
#I need to remove samples that don't have sex information for beta diversity analyses
dim(FungTableReduced)
FungTableReduced2 <- FungTableReduced[,-which(MetadataReduced$Sex=="Unknown")]
dim(FungTableReduced2)

#I should drop them from the metadata
MetadataReduced2 <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced2)

#Now I want to create a single subsample for beta diversity analyses
set.seed(1210)

SubsampleBeta<-rrarefy(t(FungTableReduced2),sample = 5000)
dim(SubsampleBeta)
rowSums(SubsampleBeta)
dim(SubsampleBeta)

#I should remove any ASVs that disappeared as I removed samples and subsampled
#Having empty ASVs will throw off the dissimilarity calculations
dim(SubsampleBeta)
ZeroList <- which(colSums(SubsampleBeta)==0)
SubsampleBeta <- SubsampleBeta[,-which(colSums(SubsampleBeta)==0)]
dim(SubsampleBeta)
#Subsample
#Dropping two samples that have extreme values in NMDS plots where they are included
#See Supplemental Table 2
SubsampleBeta <- SubsampleBeta[-c(43,47),]
dim(SubsampleBeta)
dim(MetadataReduced)
MetadataReduced2 <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced2)
MetadataReduced2 <- MetadataReduced2[-c(43,47),]
dim(MetadataReduced2)
#I am going to create a Bray-Curtis dissimilarity matrix
BrayDist<-vegdist(SubsampleBeta,method="bray",binary=FALSE)
#Here is a Jaccard distance matrix
JaccDist<- vegdist(SubsampleBeta,method="jaccard",binary=TRUE)


#Now the data are ready for analyses
#Let's start with the PERMANOVA of the Bray-Curtis dissimilarity matrix
BrayPERMANOVA<-adonis2(BrayDist~ MetadataReduced2$Age + MetadataReduced2$Sex + MetadataReduced2$Season + MetadataReduced2$Group + MetadataReduced2$Preservative,by="margin")

BrayPERMANOVA
#Season, preservative, and group are significant in the PERMANOVA

#Let's run NMDS
BrayNMDS<-metaMDS(BrayDist, k = 2,autotransform=FALSE,try=100,trymax=100)

#I want to look at the beta dispersion of samples by seasons
SeasonBrayBetaDisper <- betadisper(d=BrayDist,group=MetadataReduced2$Season)
anova(SeasonBrayBetaDisper)
TukeyHSD(SeasonBrayBetaDisper)

#Now for sex
SexBrayBetaDisper <- betadisper(d=BrayDist,group=MetadataReduced2$Sex)
anova(SexBrayBetaDisper)
#The dispersion of points is not significantly different


#Now to plot the Bray-Curtis distances

#This creates a base plot with gg_ordiplot, which gives us access to specific desireable ellipses
PlotEBase <- gg_ordiplot(ord=BrayNMDS$points, groups = MetadataReduced2$Season, kind = "sd", ellipse = TRUE)

#This plot is more flexible and interfaces with cowplot better than the gg_ordiplot
PlotE <- ggplot(data = PlotEBase$df_ord, aes(x = x, y = y, color = Group)) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  annotate("text",x=-1,y=2.47,label="italic(p) - value == 0.048", parse=TRUE, hjust =0 , size=6) +
  annotate("text",x=-1,y=2.26,label=expression( R ^"2"~ "= 0.028"), hjust =0, size=6) +
  annotate("text",x=-1,y=1.97,label="Stress = 0.158" , hjust =0, size=6) +
  geom_point(size = 3, show.legend=FALSE) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(aspect.ratio = 1) +
  xlim(c(-1.07,2.5)) +
  ylim(c(-1.07, 2.5)) +
  geom_path(data = PlotEBase$df_ellipse, show.legend=FALSE)
PlotE

PlotE2Base <- gg_ordiplot(ord=BrayNMDS$points, groups = MetadataReduced2$Sex, kind = "sd", ellipse = TRUE)

PlotE2 <- ggplot(data = PlotE2Base$df_ord, aes(x = x, y = y, color = Group)) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  annotate("text",x=-1,y=2.47,label="italic(p) - value == 0.358", parse=TRUE, hjust =0, size=6) +
  annotate("text",x=-1,y=2.26,label=expression( R ^"2"~ "= 0.020"), hjust =0, size = 6) +
  annotate("text",x=-1,y=1.97,label="Stress = 0.158" , hjust =0, size = 6) +
  geom_point(size = 3, show.legend=FALSE) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(aspect.ratio = 1) +
  xlim(c(-1.07,2.5)) +
  ylim(c(-1.07, 2.5)) +
  geom_path(data = PlotE2Base$df_ellipse, show.legend=FALSE)
PlotE2

#Let's move on to the Jaccard distance PERMANOVA
JaccPERMANOVA<-adonis2(JaccDist~ MetadataReduced2$Age + MetadataReduced2$Sex + MetadataReduced2$Season + MetadataReduced2$Group + MetadataReduced2$Preservative,by="margin")

JaccPERMANOVA
#Season, Group, and Preservative are significant

#Let's run NMDS
JaccNMDS<-metaMDS(JaccDist, k = 2,autotransform=FALSE,try=100,trymax=100)
JaccNMDS

PlotJaccABase <- gg_ordiplot(ord=JaccNMDS$points, groups = MetadataReduced2$Season, kind = "sd", ellipse = TRUE)

PlotJaccA <- ggplot(data = PlotJaccABase$df_ord, aes(x = x, y = y, color = Group)) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  scale_color_manual(values=c("#CC79A7","#009E73")) +
  annotate("text",x=0.2,y=0.45,label="italic(p) - value == 0.008", parse=TRUE, hjust =0, size = 6) +
  annotate("text",x=0.2,y=0.4,label=expression( R ^"2"~ "= 0.024"), hjust =0, size = 6) +
  annotate("text",x=0.2,y=0.35,label="Stress = 0.178" , hjust =0, size = 6) +
  geom_point(size = 3, show.legend=TRUE) +
  labs(color = "Climatic Period") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(aspect.ratio = 1) +
  geom_path(data = PlotJaccABase$df_ellipse, show.legend=FALSE)
PlotJaccA

PlotJaccBBase <- gg_ordiplot(ord=JaccNMDS$points, groups = MetadataReduced2$Sex, kind = "sd", ellipse = TRUE)

PlotJaccB <- ggplot(data = PlotJaccBBase$df_ord, aes(x = x, y = y, color = Group, shape = Group)) +
  theme_bw() +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  annotate("text",x=0.3,y=0.45,label="italic(p) - value == 0.128", parse=TRUE, hjust =0) +
  annotate("text",x=0.3,y=0.4,label=expression( R ^"2"~ "= 0.021"), hjust =0) +
  annotate("text",x=0.3,y=0.35,label="Stress = 0.178" , hjust =0) +
  geom_point(size = 3, show.legend=FALSE) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(aspect.ratio = 1) +
  geom_path(data = PlotJaccBBase$df_ellipse, show.legend=FALSE)
PlotJaccB

# Generating Figure 2----
#I exported this with a width of 800 pixels an a height of 1320 pixels
plot_grid(PlotB2 + theme(legend.position="none"),
          PlotB + theme(legend.position="none"),
          PlotC2,PlotC,PlotE2,PlotE,legendSex,legendSeason,nrow=4,ncol=2,
          rel_heights = c(1,1,1,.3),
          rel_widths = c(1,1),
          labels = c("A","B","C","D","E","F"), align = "hv")

# Loading and preparing diet data----
#Setting the working directory
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData")

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

dim(RarefiedAlphaDiversityTable)
RarefiedAlphaDiversityTable <- RarefiedAlphaDiversityTable[,-which(MetadataReduced$Sex == "Unknown")]
dim(RarefiedAlphaDiversityTable)

dim(MetadataReduced)
MetadataReduced <- MetadataReduced[-which(MetadataReduced$Sex=="Unknown"),]
dim(MetadataReduced)

#I need to match these up to the diet data
MetadataReduced$sampleid
colnames(Plant)[3:180]
#There are approximately triple the number of samples that the mycobiome has at this point
#I am going to pull out the plant information for the samples that I have in my mycobiome work

#I am going to make a table of the 60 samples that have sufficient fungal reads
#This will inform mantel tests comparing diet to fungi
Plant60 <- Plant[,which(metadata$sampleid[1]==colnames(Plant))]
SelectListPlant60 <- which(metadata$sampleid[1]==colnames(Plant))
for(i in 2:length(metadata$sampleid)){
  Plant60 <- cbind(Plant60,Plant[,which(metadata$sampleid[i]==colnames(Plant))])
  SelectListPlant60 <- c(SelectListPlant60,which(metadata$sampleid[i]==colnames(Plant)))
}

dim(Plant60)
colnames(Plant60) <- colnames(Plant)[SelectListPlant60]
#Now to drop the low read samples
Plant60 <- Plant60[,-which(colSums(FungTable)<5000)]
dim(Plant60)

#This will pull out the 52 samples that we have sufficient reads and full sex information
PlantsReduced <- Plant[,which(MetadataReduced$sampleid[1]==colnames(Plant))]
SelectListPlant <- which(MetadataReduced$sampleid[1]==colnames(Plant))
for(i in 2:length(MetadataReduced$sampleid)){
  PlantsReduced <- cbind(PlantsReduced,Plant[,which(MetadataReduced$sampleid[i]==colnames(Plant))])
  SelectListPlant <- c(SelectListPlant,which(MetadataReduced$sampleid[i]==colnames(Plant)))
}

dim(PlantsReduced)

colnames(PlantsReduced) <- colnames(Plant)[SelectListPlant]

# Calculating within sample diversity metrics of the marmoset diet----
#I am going to make the alpha diversity metrics for these samples
#Here are alpha diversity metrics for the plants
ShannonDiversityPlants <- diversity(t(PlantsReduced), index="shannon")
SimpsonDiversityPlants <- diversity(t(PlantsReduced), index="simpson")
Chao1DiversityPlants <- apply(t(PlantsReduced), MARGIN = 1,chao1, taxa.row=FALSE)
RichnessPlants <- rowSums(ifelse(t(PlantsReduced)>0,1,0))
PielouPlants <- ShannonDiversityPlants/log(RichnessPlants)

#Let's get my set of 60 samples for the invertebrates
Inverts60 <- Invert[,which(metadata$sampleid[1]==colnames(Invert))]
SelectListInvert60 <- which(metadata$sampleid[1]==colnames(Invert))
for(i in 2:length(metadata$sampleid)){
  Inverts60 <- cbind(Inverts60,Invert[,which(metadata$sampleid[i]==colnames(Invert))])
  SelectListInvert60 <- c(SelectListInvert60,which(metadata$sampleid[i]==colnames(Invert)))
}

dim(Inverts60)
colnames(Inverts60) <- colnames(Invert)[SelectListInvert60]
Inverts60
#Now to drop the low read samples
Inverts60 <- Inverts60[,-which(colSums(FungTable)<5000)]
dim(Inverts60)

#Let's get my set of 52 samples for the invertebrates
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

# Loading and preparing results of FUNGuild analysis----
#Setting a working directory to load FUNGuild analysis results
setwd("C:/Users/LabTop01/Box/GombashDocuments/MarmosetData/ITSWork/Qiime2Outputs/FUNGuild")

FunGuild<- read.table("ASVTable.taxa.guilds.txt",sep="\t",header=TRUE)
length(unique(FunGuild$guild))
length(grep("Animal Pathogen",FunGuild$guild))
length(grep("Animal Endosymbiont",FunGuild$guild))
length(grep("na", FunGuild$guild))
#Of the 2099 ASVs, only 818 have any guild information, 1281 have no information

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

#I want to drop things that are no longer in the FungTableReduced
FungusKeyReduced <- FungusKey[-which(rowSums(FungTableReduced2)==0),]

#Let's get rid of things that no longer exist
length(which(rowSums(FungTableReduced2)==0))

TaxTableReduced <- TaxTable[-which(rowSums(FungTableReduced2)==0),]
dim(TaxTableReduced)

FungTableReduced3 <- FungTableReduced2[-which(rowSums(FungTableReduced2)==0),]
dim(FungTableReduced3)

rownames(MetadataReduced) <- MetadataReduced$sampleid

# Calculating within sample diversity metrics of specific groups of fungi----
#I want to calculate within sample diversity metrics associated with resident and transient Fungi
#There are ASVs that have ecological roles that could make them resident and transient fungi
#I am going to treat fungi that could be considered both as transient fungi 
#This is in line with literature that most fungi would be transient fungi

ShannonTableBothAsPlant<- matrix(nrow=1000,ncol=52)
RichnessTableBothAsPlant<- matrix(nrow=1000,ncol=52)
ChaoTableBothAsPlant<- matrix(nrow=1000,ncol=52)
SimpsonTableBothAsPlant<- matrix(nrow=1000,ncol=52)
PielouTableBothAsPlant<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTableBothAsPlant)){
  Subsample<-rrarefy(t(FungTableReduced3),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,2]==1)]
  ShannonTableBothAsPlant[i,]<-diversity(Subsample, index="shannon")
  SimpsonTableBothAsPlant[i,]<-diversity(Subsample, index="simpson")
  ChaoTableBothAsPlant[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTableBothAsPlant[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTableBothAsPlant[i,]<- ShannonTableBothAsPlant[i,]/log(RichnessTableBothAsPlant[i,])
}

RarefiedAlphaDiversityTableBothAsPlant <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableBothAsPlant) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableBothAsPlant) <- colnames(FungTableReduced2)
RarefiedAlphaDiversityTableBothAsPlant[1,] <- colMeans(ShannonTableBothAsPlant)
RarefiedAlphaDiversityTableBothAsPlant[2,] <- colMeans(RichnessTableBothAsPlant)
RarefiedAlphaDiversityTableBothAsPlant[3,] <- colMeans(ChaoTableBothAsPlant)
RarefiedAlphaDiversityTableBothAsPlant[4,] <- colMeans(SimpsonTableBothAsPlant)
RarefiedAlphaDiversityTableBothAsPlant[5,] <- colMeans(PielouTableBothAsPlant)
RarefiedAlphaDiversityTableBothAsPlant

#Now I want to only count the ones that are just resident and not also transient

ShannonTableOnlyAniPath<- matrix(nrow=1000,ncol=52)
RichnessTableOnlyAniPath<- matrix(nrow=1000,ncol=52)
ChaoTableOnlyAniPath<- matrix(nrow=1000,ncol=52)
SimpsonTableOnlyAniPath<- matrix(nrow=1000,ncol=52)
PielouTableOnlyAniPath<- matrix(nrow=1000,ncol=52)

set.seed(50)

for (i in 1:nrow(ShannonTableOnlyAniPath)){
  Subsample<-rrarefy(t(FungTableReduced3),sample = 5000)
  Subsample<-Subsample[,which(FungusKeyReduced[,2]==0 & FungusKeyReduced[,1]==1)]
  ShannonTableOnlyAniPath[i,]<-diversity(Subsample, index="shannon")
  SimpsonTableOnlyAniPath[i,]<-diversity(Subsample, index="simpson")
  ChaoTableOnlyAniPath[i,]<- apply(Subsample, MARGIN = 1,chao1, taxa.row=FALSE)
  RichnessTableOnlyAniPath[i,]<- rowSums(ifelse(Subsample>0,1,0))
  PielouTableOnlyAniPath[i,]<- ShannonTableOnlyAniPath[i,]/log(RichnessTableOnlyAniPath[i,])
}


RarefiedAlphaDiversityTableOnlyAniPath <- matrix(nrow = 5, ncol = 52)
rownames(RarefiedAlphaDiversityTableOnlyAniPath) <- c("Shannon", "Richness", "Chao1", "Simpson", "Pielou")
colnames(RarefiedAlphaDiversityTableOnlyAniPath) <- colnames(FungTableReduced2)
RarefiedAlphaDiversityTableOnlyAniPath[1,] <- colMeans(ShannonTableOnlyAniPath)
RarefiedAlphaDiversityTableOnlyAniPath[2,] <- colMeans(RichnessTableOnlyAniPath)
RarefiedAlphaDiversityTableOnlyAniPath[3,] <- colMeans(ChaoTableOnlyAniPath)
RarefiedAlphaDiversityTableOnlyAniPath[4,] <- colMeans(SimpsonTableOnlyAniPath)
RarefiedAlphaDiversityTableOnlyAniPath[5,] <- colMeans(PielouTableOnlyAniPath)
RarefiedAlphaDiversityTableOnlyAniPath

# ANCOM-BC Analysis----
#ANCOM-BC does not require us to drop low read samples, it accounts for that already
#I am dropping samples that are missing relevant metadata
FungTableANCOMBC <- FungTable[,-which(metadata$Sex=="Unknown")]
dim(FungTableANCOMBC)
MetadataANCOMBC <- metadata[-which(metadata$Sex=="Unknown"),]
dim(MetadataANCOMBC)
rownames(MetadataANCOMBC) <- MetadataANCOMBC$sampleid
#FungTableANCOMBC <- FungTableANCOMBC[-which(rowSums(FungTableANCOMBC==0)),]

#Now let's make a phyloseq object, which is needed to run ANCOM-BC
FungPhylo <- phyloseq(otu_table(FungTableANCOMBC,taxa_are_rows = TRUE),tax_table(as.matrix(TaxTable)),sample_data(MetadataANCOMBC))

#I need to make a tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(FungPhylo)
#Need to figure out if the patch thing still exists
FungANCOM <- ancombc2(data = tse, fix_formula = "Sex + Season", p_adj_method = "BH")
FungANCOM$res

ANCOMSexResults <- data.frame(FungANCOM$res[which(FungANCOM$res$diff_SexMale==TRUE),c(1,3,6)])
ANCOMSexNames <- c("Fungus 1", "Nothophoma sp.", "Cladosporium 1","Leptosillia sp.", "Cladosporium 2", "Fungus 2", "Fungus 3")
ANCOMSexResults <- data.frame(cbind(ANCOMSexResults, ANCOMSexNames))

ANCOMSeasonResults <- data.frame(FungANCOM$res[which(FungANCOM$res$diff_SeasonWet==TRUE),c(1,4,7)])
ANCOMSeasonNames <- c("Nothophoma sp.", "Cladosporium 1", "Cladosporium 2", "Candida parapsilosis", "Penicillium sumatraense")
ANCOMSeasonResults <- data.frame(cbind(ANCOMSeasonResults,ANCOMSeasonNames))

#Preparing the data for the figure

#This stuff helps me italicize some of the names
Names1 <- ANCOMSexResults$ANCOMSexNames
Italics1 <- str_starts(Names1, "Fungus", negate=TRUE)
StyledClass1 <- ifelse(Italics1==TRUE,glue("<i>{Names1}</i>"), Names1)

ANCOMSexResults <- cbind(StyledClass1, ANCOMSexResults)

Names2 <- ANCOMSeasonResults$ANCOMSeasonNames
Italics2 <- str_starts(Names2, "Fungus", negate=TRUE)
StyledClass2 <- ifelse(Italics2==TRUE, glue("<i>{Names2}</i>"), Names2)

ANCOMSeasonResults <- cbind(StyledClass2, ANCOMSeasonResults)

#ANCOMBC Plots
PlotG <- ggplot(ANCOMSexResults, aes(y=StyledClass1, x = lfc_SexMale)) +
  theme_bw() +
  geom_point(show.legend = FALSE, size = 4,color = c("#56B4E9","#56B4E9","#E69F00","#56B4E9","#E69F00","#56B4E9","#56B4E9")) +
  theme( text = element_text(size = 17),axis.text.y = ggtext::element_markdown()) +
  geom_segment(aes(x=(lfc_SexMale + se_SexMale),y=StyledClass1,xend=(lfc_SexMale - se_SexMale),yend=StyledClass1),inherit.aes = FALSE, linewidth=1.2,color = c("#56B4E9","#56B4E9","#E69F00","#56B4E9","#E69F00","#56B4E9","#56B4E9")) +
  geom_segment(aes(x=0,xend=0,y=0,yend=7.6)) +
  xlab("Logfold Change") +
  xlim(-5.6,4.25) +
  ylab("Taxon")
PlotG

PlotH <- ggplot(ANCOMSeasonResults, aes(y=StyledClass2, x = lfc_SeasonWet)) +
  theme_bw() +
  geom_point(show.legend = FALSE, size = 4,color = c("#CC79A7","#009E73","#009E73","#009E73","#009E73")) +
  theme( text = element_text(size = 17),axis.text.y = ggtext::element_markdown()) +
  geom_segment(aes(x=(lfc_SeasonWet + se_SeasonWet),y=StyledClass2,xend=(lfc_SeasonWet - se_SeasonWet),yend=StyledClass2),inherit.aes = FALSE, linewidth=1.2,color = c("#CC79A7","#009E73","#009E73","#009E73","#009E73")) +
  geom_segment(aes(x=0,xend=0,y=0,yend=5.6)) +
  xlab("Logfold Change") +
  xlim(-5.6,4.25) +
  ylab("Taxon")
PlotH
# Generating Figure 3----
#Export at width 800 pixels and height 800 pixels
plot_grid(PlotG + theme(legend.position="none"),
          legendSex,
          PlotH + theme(legend.position="none"),
          legendSeason,nrow=4,ncol=1,
          rel_heights = c(1,.3,1,.3),
          labels = c("A","","B",""), align = "hv")


# Correlation analyses between diet and fungi----
PlotObject<-data.frame(t(rbind(
  RarefiedAlphaDiversityTable,
  RarefiedAlphaDiversityTableOnlyAniPath,
  RarefiedAlphaDiversityTableBothAsPlant,
  SimpsonDiversityInverts,
  SimpsonDiversityPlants,
  ShannonDiversityInverts,
  ShannonDiversityPlants,
  RichnessPlants,
  RichnessInverts,
  PielouPlants,
  PielouInverts,
  Chao1DiversityInverts,
  Chao1DiversityPlants)))

#Running correlation tests between the overall mycobiome and the two diet components
MycoRichVsInsectRich <- cor.test(x = PlotObject[,2], y = PlotObject[,21], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoRichVsInsectRich

MycoRichVsPlantRich <- cor.test(x = PlotObject[,2], y = PlotObject[,20], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoRichVsPlantRich

MycoShannonVsInsectShannon <- cor.test(x = PlotObject[,1], y = PlotObject[,18], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoShannonVsInsectShannon

MycoShannonVsPlantShannon <- cor.test(x = PlotObject[,1], y = PlotObject[,19], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoShannonVsPlantShannon

MycoChaoVsInsectChao <- cor.test(x = PlotObject[,3], y = PlotObject[,24], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoChaoVsInsectChao

MycoChaoVsPlantChao <- cor.test(x = PlotObject[,3], y = PlotObject[,25], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoChaoVsPlantChao

MycoSimpsonVsInsectSimpson <- cor.test(x = PlotObject[,4], y = PlotObject[,16], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoSimpsonVsInsectSimpson

MycoSimpsonVsPlantSimpson <- cor.test(x = PlotObject[,4], y = PlotObject[,17], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoSimpsonVsPlantSimpson

MycoPielouVsInsectPielou <- cor.test(x = PlotObject[,5], y = PlotObject[,23], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoPielouVsInsectPielou

MycoPielouVsPlantPielou <- cor.test(x = PlotObject[,5], y = PlotObject[,22], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
MycoPielouVsPlantPielou

#Now for the likely resident taxa
ResidentRichVsInsectRich <- cor.test(x = PlotObject[,7], y = PlotObject[,21], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentRichVsInsectRich

ResidentRichVsPlantRich <- cor.test(x = PlotObject[,7], y = PlotObject[,20], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentRichVsPlantRich

ResidentShannonVsInsectShannon <- cor.test(x = PlotObject[,6], y = PlotObject[,18], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentShannonVsInsectShannon

ResidentShannonVsPlantShannon <- cor.test(x = PlotObject[,6], y = PlotObject[,19], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentShannonVsPlantShannon

ResidentChaoVsInsectChao <- cor.test(x = PlotObject[,8], y = PlotObject[,24], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentChaoVsInsectChao

ResidentChaoVsPlantChao <- cor.test(x = PlotObject[,8], y = PlotObject[,25], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentChaoVsPlantChao

ResidentSimpsonVsInsectSimpson <- cor.test(x = PlotObject[,9], y = PlotObject[,16], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentSimpsonVsInsectSimpson

ResidentSimpsonVsPlantSimpson <- cor.test(x = PlotObject[,9], y = PlotObject[,17], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentSimpsonVsPlantSimpson

ResidentPielouVsInsectPielou <- cor.test(x = PlotObject[,10], y = PlotObject[,23], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentPielouVsInsectPielou

ResidentPielouVsPlantPielou <- cor.test(x = PlotObject[,10], y = PlotObject[,22], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
ResidentPielouVsPlantPielou

#Now for the Likely Transient fungi
TransientRichVsInsectRich <- cor.test(x = PlotObject[,12], y = PlotObject[,21], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientRichVsInsectRich

TransientRichVsPlantRich <- cor.test(x = PlotObject[,12], y = PlotObject[,20], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientRichVsPlantRich

TransientShannonVsInsectShannon <- cor.test(x = PlotObject[,11], y = PlotObject[,18], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientShannonVsInsectShannon

TransientShannonVsPlantShannon <- cor.test(x = PlotObject[,11], y = PlotObject[,19], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientShannonVsPlantShannon

TransientChaoVsInsectChao <- cor.test(x = PlotObject[,13], y = PlotObject[,24], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientChaoVsInsectChao

TransientChaoVsPlantChao <- cor.test(x = PlotObject[,13], y = PlotObject[,25], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientChaoVsPlantChao

TransientSimpsonVsInsectSimpson <- cor.test(x = PlotObject[,14], y = PlotObject[,16], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientSimpsonVsInsectSimpson

TransientSimpsonVsPlantSimpson <- cor.test(x = PlotObject[,14], y = PlotObject[,17], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientSimpsonVsPlantSimpson

TransientPielouVsInsectPielou <- cor.test(x = PlotObject[,15], y = PlotObject[,23], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientPielouVsInsectPielou

TransientPielouVsPlantPielou <- cor.test(x = PlotObject[,15], y = PlotObject[,22], alternative = c("two.sided"), method = "spearman", exact = FALSE,data=PlotObject)
TransientPielouVsPlantPielou

#Now I want to make plots of these correlations

CorPlot1 <- ggplot(PlotObject, aes(y=Richness, x=RichnessInverts)) +
  geom_point(size=2, shape=19) +
  geom_smooth(method="lm",col="black") +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  xlab("Arthropod Richness") +
  ylab("Mycobiome Richness")
CorPlot1

CorPlot2 <- ggplot(PlotObject, aes(y=Richness.2, x=RichnessInverts)) +
  geom_point(size=2, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  xlab("Arthropod Richness") +
  ylab("Likely Transient Fungi Richness")

CorPlot3 <- ggplot(PlotObject, aes(y=Simpson.1, x=SimpsonDiversityInverts)) +
  geom_point(size=2, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  xlab("Arthropod Simpson Diversity") +
  ylab("Likely Resident Fungi Simpson Div.")

CorPlotA <- ggplot(PlotObject, aes(y=Richness, x=RichnessInverts)) +
  geom_point(size=3, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  geom_smooth(method="lm",col="black") +
  xlab("") +
  annotate("text",x=35,y=135,label="italic(p) - value == 0.039", parse=TRUE, hjust =0, size = 6) +
  annotate("text",x=35,y=125,label= expression( rho ~ "= 0.287" ), hjust =0, size = 6) +
  ylab("Mycobiome Richness")

CorPlotB <- ggplot(PlotObject, aes(y=Richness.2, x=RichnessInverts)) +
  geom_point(size=3, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  xlab("") +
  geom_smooth(method="lm",col="black") +
  annotate("text",x=35,y=29,label="italic(p) - value == 0.035", parse=TRUE, hjust =0, size = 6) +
  annotate("text",x=35,y=27,label= expression( rho ~ "= 0.293" ), hjust =0, size = 6) +
  ylab("Transient Fungi Richness")

CorPlotC <- ggplot(PlotObject, aes(y=Richness.1, x=RichnessInverts)) +
  geom_point(size=3, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 17)) +
  xlab("Arthropod Richness") +
  geom_smooth(method="lm",col="black") +
  annotate("text",x=35,y=5.8,label=expression(italic('p') * " - value" * " = 0.150"), hjust =0, size = 6) +
  annotate("text",x=35,y=5.4,label= expression( rho ~ "= 0.203" ), hjust =0, size = 6) +
  ylab("Resident Fungi Richness")

# Generating Figure 4----
#I exported this with a width of 400 pixels and a height of 1200 pixels
plot_grid(CorPlotA, CorPlotB, CorPlotC,nrow=3,ncol=1,labels=c("A","B","C"), rel_widths = c(1), rel_heights = c(1,1,1))

CorPlotDTalk <- ggplot(PlotObject, aes(y=Simpson.1, x=SimpsonDiversityInverts)) +
  geom_point(size=2, shape=19) +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  xlab("Arthropod Simpson Div.") +
  geom_smooth(method="lm",col="black") +
  annotate("text",x=0.7,y=0.85,label="italic(p) - value == 0.032", parse=TRUE, hjust =0, size=13) +
  annotate("text",x=0.7,y=0.75,label= expression( rho ~ "= -0.298" ), hjust =0, size=13) +
  ylab("Likely Resident Fungi Simpson Div.")
CorPlotDTalk

# Fixed and Mixed Effects Linear Regressions Testing Correlations between sex/season and transient fungi----
#Do the within sample diversity metrics of the transient fungi vary by sex?
#Now I need to see what predicts Plant Associated Alpha Diversity when I count the ASVs that are both animal and plant associated as plant associated
ModelDataTransient <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableBothAsPlant)))
#I want to make sure the factors are coded as factors
ModelDataTransient$Sex <-as.factor(ModelDataTransient$Sex)
ModelDataTransient$Group <-as.factor(ModelDataTransient$Group)
ModelDataTransient$Season <-as.factor(ModelDataTransient$Season)
ModelDataTransient$Reproductive <-as.factor(ModelDataTransient$Reproductive)
ModelDataTransient$Individual <-as.factor(ModelDataTransient$Individual)
ModelDataTransient$Preservative <-as.factor(ModelDataTransient$Preservative)
ModelDataTransient$Age <-as.factor(ModelDataTransient$Age)
ModelDataTransient$Time2 <-as.factor(ModelDataTransient$Time2)
#Now I want to put in NAs when appropriate
ModelDataTransient$Reproductive[which(ModelDataTransient$Reproductive=="Unknown")]<-NA
ModelDataTransient$Age[which(ModelDataTransient$Age=="Unknown")]<-NA

#Shannon
ShannonModelFullTransient <- lm(Shannon ~ 
                         Season +
                         Sex +
                         Group +
                         Preservative +
                         Age,
                       data = ModelDataTransient)
summary(ShannonModelFullTransient)

#Richness
RichnessModelFullTransient <- lmer(Richness ~ 
                            Season +
                            Sex +
                            (1|Group) +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelDataTransient, REML=FALSE)
summary(RichnessModelFullTransient)

RichnessModelNullSeasonTransient <- lmer(Richness ~ 
                                  #Season +
                                  Sex +
                                  (1|Group) +
                                  Preservative +
                                  Age +
                                  Season:Time2,
                                data = ModelDataTransient, REML=FALSE)
summary(RichnessModelNullSeasonTransient)
RichnessModelNullSexTransient <- lmer(Richness ~ 
                               Season +
                               #Sex +
                               (1|Group) +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelDataTransient, REML=FALSE)
summary(RichnessModelNullSexTransient)
#Now for significance tests
anova(RichnessModelNullSeasonTransient, RichnessModelFullTransient, test="Chisq")
anova(RichnessModelNullSexTransient, RichnessModelFullTransient, test="Chisq")
#Sex is significant, season is not
PlotATransient <- ggplot(ModelData, aes(x = Sex, y = Richness, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  ylab("Richness")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotATransient
#A Wald's test to see if the sexes significantly differ
Anova(RichnessModelFullTransient)

#Chao1
Chao1ModelFullTransient <- lmer(Chao1 ~ 
                         Season +
                         Sex +
                         (1|Group) +
                         Preservative +
                         Age +
                         Season:Time2,
                       data = ModelDataTransient, REML=FALSE)
summary(Chao1ModelFullTransient)

Chao1ModelNullSeasonTransient <- lmer(Chao1 ~ 
                               #Season +
                               Sex +
                               (1|Group) +
                               Preservative +
                               Age +
                               Season:Time2,
                             data = ModelDataTransient, REML=FALSE)
summary(Chao1ModelNullSeasonTransient)
Chao1ModelNullSexTransient <- lmer(Chao1 ~ 
                            Season +
                            #Sex +
                            (1|Group) +
                            Preservative +
                            Age +
                            Season:Time2,
                          data = ModelDataTransient, REML=FALSE)
summary(Chao1ModelNullSexTransient)
#Now for significance tests
anova(Chao1ModelNullSeasonTransient, Chao1ModelFullTransient, test="Chisq")
anova(Chao1ModelNullSexTransient, Chao1ModelFullTransient, test="Chisq")
#A Wald's test to see if the sexes significantly differ
Anova(Chao1ModelFullTransient)

#Sex is significant, season is not
PlotX <- ggplot(ModelData, aes(x = Sex, y = Chao1, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  ylab("Chao1")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotX

#Simpson
SimpsonModelFullTransient <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         Preservative +
                         Age,
                       data = ModelDataTransient)
summary(SimpsonModelFullTransient)

#Pielou
#Note, some samples had 0 Simpson and have NaN for Pielou, those samples were removed
PielouModelFullTransient <- lm(Pielou ~ 
                        Season +
                        Sex +
                        Group +
                        Preservative +
                        Age,
                      data = ModelDataTransient[-which(ModelDataTransient$Pielou=="NaN"),])
summary(PielouModelFullTransient)

# Fixed and Mixed Effects Linear Regressions Testing Correlations between sex/season and resident fungi----
#Now I need to run models for the resident fungi
ModelDataResident <- data.frame(cbind(MetadataReduced,t(RarefiedAlphaDiversityTableOnlyAniPath)))
#I want to make sure the factors are coded as factors
ModelDataResident$Sex <-as.factor(ModelDataResident$Sex)
ModelDataResident$Group <-as.factor(ModelDataResident$Group)
ModelDataResident$Season <-as.factor(ModelDataResident$Season)
ModelDataResident$Reproductive <-as.factor(ModelDataResident$Reproductive)
ModelDataResident$Individual <-as.factor(ModelDataResident$Individual)
ModelDataResident$Preservative <-as.factor(ModelDataResident$Preservative)
ModelDataResident$Age <-as.factor(ModelDataResident$Age)
ModelDataResident$Time2 <-as.factor(ModelDataResident$Time2)
#Now I want to put in NAs when appropriate
ModelDataResident$Reproductive[which(ModelDataResident$Reproductive=="Unknown")]<-NA
ModelDataResident$Age[which(ModelDataResident$Age=="Unknown")]<-NA


ShannonModelFullResident <- lm(Shannon ~ 
                         Season +
                         Sex +
                         Group +
                         Preservative +
                         Age,
                       data = ModelDataResident)
summary(ShannonModelFullResident)
#The model itself is not significant (p-value = 0.4309)
#Sex is the only significant parameter, with male being different from female

PlotAResident <- ggplot(ModelData, aes(x = Sex, y = Shannon, color = Sex)) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme_bw() +
  ylab("Shannon diversity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend=FALSE)
PlotAResident


RichnessModelFullResident <- lm(Richness ~ 
                          Season +
                          Sex +
                          Group +
                          Preservative +
                          Age,
                        data = ModelDataResident)
summary(RichnessModelFullResident)
#The model itself is not significant (p-value = 0.882)
#No parameter is significant



Chao1ModelFullResident <- lm(Chao1 ~ 
                       Season +
                       Sex +
                       Group +
                       Preservative +
                       Age,
                     data = ModelDataResident)
summary(Chao1ModelFullResident)
#The model itself is not significant (p-value = 0.886)
#No parameter is significant

SimpsonModelFullResident <- lm(Simpson ~ 
                         Season +
                         Sex +
                         Group +
                         Preservative +
                         Age,
                       data = ModelDataResident)
summary(SimpsonModelFullResident)
#The model itself is not significant (p-value = 0.6564)
#No parameter is significant

#Note that some samples have a Pielou of NaN, and those were removed
PielouModelFullResident <- lm(Pielou ~ 
                        Season +
                        Sex +
                        Group +
                        Preservative +
                        Age,
                      data = ModelDataResident[-which(ModelDataResident$Pielou=="NaN"),])
summary(PielouModelFullResident)
#The model itself is not significant (p-value = 0.5152)
#The Road and F group groups are different from Abocarpa or whatever it is called


# Mantel tests----
#Let's make a subsamle of 60 samples for the Mantel tests
set.seed(1210)

SubsampleMantel<-rrarefy(t(FungTableReduced),sample = 5000)
dim(SubsampleMantel)
rowSums(SubsampleMantel)
dim(SubsampleMantel)
#I need to make another fungus key that matches the dimensions of the SubSampleMantel
FungusKeyReduced2 <- FungusKey[-which(colSums(SubsampleMantel)==0),]
dim(FungusKeyReduced2)
#Let's drop the fungi that don't appear
dim(SubsampleMantel)
SubsampleMantel <- SubsampleMantel[,-which(colSums(SubsampleMantel)==0)]
dim(SubsampleMantel)

#I should create a resident fungi subsample
APSubsample <- SubsampleMantel[,which(FungusKeyReduced2[,1]==1 & FungusKeyReduced2[,2]==0)]
dim(APSubsample)
#Now a transient fungi subsample
PASubsample <- SubsampleMantel[,which(FungusKeyReduced2[,2]==1)]
dim(PASubsample)

#We are loading this package here because a dependent package interfers with ANCOM-BC
library(ecodist)
#Does the overall mycobiome vary with the invertebrate portion of the diet?
ecodist::mantel(formula=(vegdist((SubsampleMantel), method="bray", binary=FALSE)) ~(vegdist(t(Inverts60), method="bray", binary=FALSE)), nperm=100000)

#Does the overall mycobiome vary with the plant portion of the diet?
ecodist::mantel(formula=(vegdist((SubsampleMantel), method="bray", binary=FALSE)) ~(vegdist(t(Plant60), method="bray", binary=FALSE)), nperm=100000)

#Some samples don't have transient taxa, so those samples must be dropped before the analysis
#The Mantel test needs the same samples in both matrices
PASubSample2 <- PASubsample[-which(rowSums(PASubsample)==0),]
Inverts60PAReduced <- Inverts60[,-which(rowSums(PASubsample)==0)]
Plant60PAReduced <- Plant60[,-which(rowSums(PASubsample)==0)]

#Does the transient mycobiome vary with the invertebrate portion of the diet?
ecodist::mantel(formula=(vegdist((PASubSample2), method="bray", binary=FALSE)) ~(vegdist(t(Inverts60PAReduced), method="bray", binary=FALSE)), nperm=100000)

#Does the transient mycobiome vary with the plant portion of the diet?
ecodist::mantel(formula=(vegdist((PASubSample2), method="bray", binary=FALSE)) ~(vegdist(t(Plant60PAReduced), method="bray", binary=FALSE)), nperm=100000)

#Let's drop samples that lack resident fungi
APSubsample2 <- APSubsample[-which(rowSums(APSubsample)==0),]
Inverts60APReduced <- Inverts60[,-which(rowSums(APSubsample)==0)]
Plant60APReduced <- Plant60[,-which(rowSums(APSubsample)==0)]

#Does the resident mycobiome vary with the invertebrate portion of the diet?
ecodist::mantel(formula=(vegdist((APSubsample2), method="bray", binary=FALSE)) ~(vegdist(t(Inverts60APReduced), method="bray", binary=FALSE)), nperm=100000)

#Does the resident mycobiome vary with the plant portion of the diet?
ecodist::mantel(formula=(vegdist((APSubsample2), method="bray", binary=FALSE)) ~(vegdist(t(Plant60APReduced), method="bray", binary=FALSE)), nperm=100000)

# CCREPE Analysis----
#I am running CCREPE
#I need relative abundance of Fungi, Plants, and Invertebrates for ccrepe
RelPlants <- matrix(nrow=dim(PlantsReduced)[1],ncol=dim(PlantsReduced)[2])
for(i in 1:dim(PlantsReduced)[1]){
  for(j in 1:dim(PlantsReduced)[2]){
    RelPlants[i,j] <- PlantsReduced[i,j]/colSums(PlantsReduced)[j]
  }
}
colSums(RelPlants)


#I need relative abundance of Fungi, Plants, and Invertebrates for CCREPE
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

RelFung <- matrix(nrow=dim(FungTableReduced2)[1],ncol=dim(FungTableReduced2)[2])
for(i in 1:dim(FungTableReduced2)[1]){
  for(j in 1:dim(FungTableReduced2)[2]){
    RelFung[i,j] <- FungTableReduced2[i,j]/colSums(FungTableReduced2)[j]
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
library(ccrepe)
set.seed(121023)
FungInvCCREPE <- ccrepe(x = t(RelFung), y = t(RelInverts),verbose=TRUE, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"),min.subj = 26, make.output.table=FALSE)

FungInvCCREPE_pvalues <- FungInvCCREPE$p.values
FungInvCCREPE_qvalues <- FungInvCCREPE$q.values
FungInvCCREPE_zstat <- FungInvCCREPE$z.stat
FungInvCCREPE_simscore <- FungInvCCREPE$sim.score


FungInvCCREPE_pvalues_NO_NA <- FungInvCCREPE_pvalues[-which(rowSums(FungInvCCREPE_pvalues, na.rm=TRUE)==0),-which(colSums(FungInvCCREPE_pvalues, na.rm=TRUE)==0)]
FungInvCCREPE_qvalues_NO_NA <- FungInvCCREPE_qvalues[-which(rowSums(FungInvCCREPE_qvalues, na.rm=TRUE)==0),-which(colSums(FungInvCCREPE_qvalues, na.rm=TRUE)==0)]
FungInvCCREPE_zstat_NO_NA <- FungInvCCREPE_zstat[-which(rowSums(FungInvCCREPE_zstat, na.rm=TRUE)==0),-which(colSums(FungInvCCREPE_zstat, na.rm=TRUE)==0)]
FungInvCCREPE_simscore_NO_NA <- FungInvCCREPE_simscore[-which(rowSums(FungInvCCREPE_simscore, na.rm=TRUE)==0),-which(colSums(FungInvCCREPE_simscore, na.rm=TRUE)==0)]


#Now for a plant CCREPE
set.seed(11623)
FungPlantCCREPE <- ccrepe(x = t(RelFung), y = t(RelPlants), verbose=TRUE, sim.score = cor, sim.score.args = list(method="spearman", use="complete.obs"),min.subj = 26, make.output.table=FALSE)


FungPlantCCREPE_pvalues <- FungPlantCCREPE$p.values
FungPlantCCREPE_qvalues <- FungPlantCCREPE$q.values
FungPlantCCREPE_zstat <- FungPlantCCREPE$z.stat
FungPlantCCREPE_simscore <- FungPlantCCREPE$sim.score


FungPlantCCREPE_pvalues_NO_NA <- FungPlantCCREPE_pvalues[-which(rowSums(FungPlantCCREPE_pvalues, na.rm=TRUE)==0),-which(colSums(FungPlantCCREPE_pvalues, na.rm=TRUE)==0)]
FungPlantCCREPE_qvalues_NO_NA <- FungPlantCCREPE_qvalues[-which(rowSums(FungPlantCCREPE_qvalues, na.rm=TRUE)==0),-which(colSums(FungPlantCCREPE_qvalues, na.rm=TRUE)==0)]
FungPlantCCREPE_zstat_NO_NA <- FungPlantCCREPE_zstat[-which(rowSums(FungPlantCCREPE_zstat, na.rm=TRUE)==0),-which(colSums(FungPlantCCREPE_zstat, na.rm=TRUE)==0)]
FungPlantCCREPE_simscore_NO_NA <- FungPlantCCREPE_simscore[-which(rowSums(FungPlantCCREPE_simscore, na.rm=TRUE)==0),-which(colSums(FungPlantCCREPE_simscore, na.rm=TRUE)==0)]
