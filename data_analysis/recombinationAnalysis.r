

# Packages
library(svglite)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)

recombinationComparisonData <- read.csv("./data/out_Recombination_LIFPNG_1750177457.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, RecombinationRate, survivalRate, loadSurvivalRate,count) %>% 
  summarise()



# plot(survivalRate+loadSurvivalRate~RecombinationRate,ylab="Surivival Rate",data=filter(recombinationComparisonData, Xlinked==0),ylim=c(0,1),log="x")
# par(new=TRUE)
# plot(loadSurvivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 19,ylim=c(0,1),log="x")
# par(new=TRUE)
# plot(survivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 2,ylim=c(0,1),log="x")
# abline(v=0.00101, col="blue")
# abline(v=1.0e-5, col="red")

# plot(survivalRate+loadedSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),ylim=c(0,1),log="x")
# par(new=TRUE)
# plot(loadSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 19,ylim=c(0,1),log="x")
# par(new=TRUE)
# plot(survivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 2,ylim=c(0,1),log="x")



# create data
autosomalRecombinationComparisonData<- filter(recombinationComparisonData, Xlinked==0, Sterile==0, PostCompetitionMutationTiming==0,GrowthRate==3.0)
RecombinationRate <- rep(autosomalRecombinationComparisonData$RecombinationRate,times=3)
outcomeProportion <- c(1-autosomalRecombinationComparisonData$survivalRate - autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$survivalRate)
outcomeGroup <- rep(c("Extinction","Pseudo-overdominance survival","Carrying capacity reached"),each=length(autosomalRecombinationComparisonData$RecombinationRate))
data <- data.frame(RecombinationRate, outcomeProportion, outcomeGroup)

# reverse stack order
data$outcomeGroup <- factor(data$outcomeGroup, c("Carrying capacity reached","Pseudo-overdominance survival","Extinction"))

outcomeProportion[autosomalRecombinationComparisonData$RecombinationRate==0]
# stacked area chart

svglite("figures/figure_4_2.svg", width = 8.2, height = 5.5)
ggplot(data,log="x", aes(x=RecombinationRate, y=outcomeProportion, fill=outcomeGroup)) + 
    geom_area() +
    scale_x_continuous(trans='log10',limits = c(min(RecombinationRate[RecombinationRate!=0]), max(RecombinationRate)),
    breaks = trans_breaks("log10", function(x) 10^(x)),expand = c(0,0)) +
    xlab("Recombination rate (M/bp)") + 
    ylab("Outcome proportion") +
    guides(fill=guide_legend(title="Outcome group")) +
    scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
    scale_y_continuous(expand = c(0, 0)) +
    annotation_logticks(sides="b") #+
    # geom_vline(xintercept=2.1e-8, color="yellow")
dev.off()

# install.packages("pracma")

require(pracma)

# Looking at where unlinked extinction is lower than linked extinction

recombinationParameterSpaceData <- read.csv("./data/out_LociSpaceComparison_LIFPNG_1750176677.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount, survivalRate, loadSurvivalRate,count) %>% 
  summarise()

# create data
slices = list()
spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.083,GrowthRate==4,MutationCount==100)
slices[[1]] = spaceSlice
spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.099,GrowthRate==4,MutationCount==70)
slices[[2]] = spaceSlice
spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.131,GrowthRate==4,MutationCount==40)
slices[[3]] = spaceSlice
# spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.2,GrowthRate==3,MutationCount==25)
spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.259,GrowthRate==4,MutationCount==10)
slices[[4]] = spaceSlice
# slices[[5]] = spaceSlice
for (spaceSlice in slices) {
  freq = spaceSlice$MutationFrequency[1]
  lociCount = spaceSlice$MutationCount[1]

  linkedData <- filter(spaceSlice,ChromosomeCount==1)
  individuals <- rep(linkedData$Individuals,times=3)
  outcomeProportion <- c(1-linkedData$survivalRate - linkedData$loadSurvivalRate,linkedData$loadSurvivalRate,linkedData$survivalRate)
  outcomeGroup <- rep(c("Extinction","Loaded survival","Survival"),each=length(linkedData$MutationCount))
  data <- data.frame(individuals, outcomeProportion, outcomeGroup)
  print("linked")
  df <- filter(data,outcomeGroup=="Extinction")
  print(trapz(df$individuals,df$outcomeProportion)/0.9) # 0.9 is the area of the space
  df <- filter(data,outcomeGroup=="Loaded survival")
  print(trapz(df$individuals,df$outcomeProportion)/0.9)
  df <- filter(data,outcomeGroup=="Survival")
  print(trapz(df$individuals,df$outcomeProportion)/0.9)


  # print(ggplot(data, aes(x=individuals, y=outcomeProportion, fill=outcomeGroup)) + 
  #     # ggtitle(sprintf("Freq %f , loci %d",freq,lociCount)) +
      
  #     theme_void() +
  #     theme(panel.background = element_blank()) +
  #     scale_x_continuous(breaks = round(seq(min(individuals), max(individuals), by = 15),1)) +
  #     xlab("Individuals") + theme(legend.position = "none") + 
  #     ylab("Outcome proportion") +
  #     geom_area() + guides(fill=guide_legend(title="Outcome group")) +
  #     annotate("text", x=20, y=0.875, label= round(AUC,2))
  #     )

  unlinkedData <- filter(spaceSlice,ChromosomeCount!=1)
  individuals <- rep(unlinkedData$Individuals,times=3)
  outcomeProportion <- c(1-unlinkedData$survivalRate - unlinkedData$loadSurvivalRate,unlinkedData$loadSurvivalRate,unlinkedData$survivalRate)
  outcomeGroup <- rep(c("Extinction","Loaded survival","Survival"),each=length(unlinkedData$MutationCount))
  data <- data.frame(individuals, outcomeProportion, outcomeGroup)
  print("unlinked")
  df <- filter(data,outcomeGroup=="Extinction")
  print(trapz(df$individuals,df$outcomeProportion)/0.9)
  df <- filter(data,outcomeGroup=="Loaded survival")
  print(trapz(df$individuals,df$outcomeProportion)/0.9)
  df <- filter(data,outcomeGroup=="Survival")
  print(trapz(df$individuals,df$outcomeProportion)/0.9)
  # print(ggplot(data, aes(x=individuals, y=outcomeProportion, fill=outcomeGroup)) + 
  #     # ggtitle(sprintf("Freq %f , loci %d",freq,lociCount)) +
  #     theme_void() +
  #     theme(panel.background = element_blank()) +
  #     scale_x_continuous(breaks = round(seq(min(individuals), max(individuals), by = 15),1)) +
  #     theme(legend.position = "none") + 
  #     xlab("Individuals") + 
  #     ylab("Outcome proportion") +
  #     guides(fill=guide_legend(title="Outcome group")) +
  #     geom_area() + 
  #     annotate("text", x=20, y=0.875, label= round(AUC,2))
  #     )
}

# install.packages("tidyverse")

library(tidyverse)
spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.083,MutationCount==100)
spaceSlices <- rbind(spaceSlice,filter(recombinationParameterSpaceData,MutationFrequency==0.099,MutationCount==70))
spaceSlices <- rbind(spaceSlices,filter(recombinationParameterSpaceData,MutationFrequency==0.131,MutationCount==40))
spaceSlices <- rbind(spaceSlices,filter(recombinationParameterSpaceData,MutationFrequency==0.259,MutationCount==10))

# spaceSlice <- filter(recombinationParameterSpaceData,MutationFrequency==0.1,MutationCount==100)
# spaceSlices <- rbind(spaceSlice,filter(recombinationParameterSpaceData,MutationFrequency==0.1,MutationCount==70))
# spaceSlices <- rbind(spaceSlices,filter(recombinationParameterSpaceData,MutationFrequency==0.15,MutationCount==40))
# spaceSlices <- rbind(spaceSlices,filter(recombinationParameterSpaceData,MutationFrequency==0.4,MutationCount==10))


slice <- filter(spaceSlices,GrowthRate==3.5)
individuals <- rep(slice$Individuals,times=3)
loci <- factor(rep(slice$MutationCount,times=3))
freq <- rep(slice$MutationFrequency,times=3)
load <- rep(paste0(slice$MutationFrequency," frequency at ",slice$MutationCount, " loci") ,times=3)
linked <- ifelse(rep(slice$ChromosomeCount,times=3)==1,"Linked","Unlinked")
outcomeProportion <- c(1-slice$survivalRate - slice$loadSurvivalRate,slice$loadSurvivalRate,slice$survivalRate)
outcomeGroup <- rep(c("Extinction","Loaded survival","Survival"),each=length(slice$Individuals))
data <- data.frame(individuals, outcomeProportion, outcomeGroup, linked, load, loci, freq)
# df <- filter(data,outcomeGroup=="Extinction")
# AUC <- trapz(df$individuals,df$outcomeProportion)

data$outcomeGroup <- factor(data$outcomeGroup, sort(unique(data$outcomeGroup), decreasing = T))
# data

      
# ggplot(data, aes(x=individuals, y=outcomeProportion, fill=outcomeGroup)) + 
#       # ggtitle(sprintf("Freq %f , loci %d",freq,lociCount)) +
#       scale_x_continuous(breaks = round(seq(min(individuals), max(individuals), by = 15),1)) +
#       xlab("Individuals") + 
#       theme(legend.position = "none") + 
#       ylab("Outcome proportion") +
#       geom_area() + guides(fill=guide_legend(title="Outcome group"))+
#       facet_grid(vars(load), vars(linked)) 

svglite("figures/figure_S42.svg", width = 4.1, height = 5.5)

ggplot(data, aes(x=individuals, y=outcomeProportion, fill=outcomeGroup)) + 
      # ggtitle(sprintf("Freq %f , loci %d",freq,lociCount)) +
      scale_x_continuous(breaks = round(seq(min(individuals), max(individuals), by = 15),1)) +
      xlab("Individuals") + 
      theme(panel.background = element_blank()) +
      theme(legend.position = "none") + 
      ylab("Outcome proportion") +
      scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
      geom_area() + 
      # geom_text(data = ann_text,label = "Text") +
      guides(fill=guide_legend(title="Outcome group"))+
      facet_grid(vars(fct_rev(loci),freq), vars(linked))
dev.off()