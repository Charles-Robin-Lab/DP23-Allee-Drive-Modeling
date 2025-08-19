

# Packages
library(svglite)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)

recombinationComparisonData <- read.csv("./data/out_Recombination_LIFPNG_1752543361.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, RecombinationRate, survivalRate, loadSurvivalRate,count) %>% 
  summarise()

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

svglite("figures/figure_6.svg", width = 8.2, height = 5.5)
ggplot(data,log="x", aes(x=RecombinationRate, y=outcomeProportion, fill=outcomeGroup)) + 
    geom_area() +
    scale_x_continuous(trans='log10',limits = c(min(RecombinationRate[RecombinationRate!=0]), max(RecombinationRate)),
    breaks = trans_breaks("log10", function(x) 10^(x)),expand = c(0,0)) +
    xlab("Recombination rate (M/bp)") + 
    ylab("Outcome proportion") +
    guides(fill=guide_legend(title="Outcome group")) +
    scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
    scale_y_continuous(expand = c(0, 0)) +
    annotation_logticks(sides="b") 
dev.off()

# install.packages("pracma")
require(pracma)

# Looking at where unlinked extinction is lower than linked extinction

recombinationParameterSpaceData <- read.csv("./data/out_LociSpaceComparison_LIFPNG_1752543377.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount, survivalRate, loadSurvivalRate,count) %>% 
  summarise()

# create data
recombinationParameterSpaceDataWithGrowthrate <- filter(recombinationParameterSpaceData,GrowthRate==3)
loads <- data.frame(
  mf = rep(c(0.083,0.099,0.131,0.259),each=2),
  mc = rep(c(100,70,40,10),each=2),
  linked = rep(c(TRUE,FALSE),times=4)
)
same_Le_datas<-apply(loads, 1, function(load) {
  data<-filter(recombinationParameterSpaceDataWithGrowthrate,
    MutationFrequency==load[["mf"]],
    MutationCount==load[["mc"]],
    (ChromosomeCount==1)==load[["linked"]])
  
  data
})
same_Le_data <- do.call(rbind, same_Le_datas)

library(tidyverse)

outcomeGroup <- rep(c("Extinction","Loaded survival","Survival"),each=length(same_Le_data$Individuals))
outcomeProportion <- c(1-same_Le_data$survivalRate - same_Le_data$loadSurvivalRate,same_Le_data$loadSurvivalRate,same_Le_data$survivalRate)
linked <- ifelse(rep(same_Le_data$ChromosomeCount,times=3)==1,"Linked","Unlinked")
individuals <- rep(same_Le_data$Individuals,times=3)
loci <- factor(rep(same_Le_data$MutationCount,times=3))
freq <- rep(same_Le_data$MutationFrequency,times=3)
data <- data.frame(individuals, outcomeProportion, outcomeGroup, linked, loci, freq)
# otherwise gg plot uses a weird ordering
data$outcomeGroup <- factor(data$outcomeGroup, sort(unique(data$outcomeGroup), decreasing = TRUE))



apply(loads,1,function(load) {
  print(sprintf("freq=%f lociCount=%f linked=%s",load[["mf"]],load[["mc"]],load[["linked"]]))
  for (group in c("Extinction","Loaded survival","Survival")) {
    df <- filter(data,
    outcomeGroup==group,
    freq==load[["mf"]],
    loci==load[["mc"]],
    (linked=="Linked")==load[["linked"]])
    print(trapz(df$individuals,df$outcomeProportion)/0.9) # 90 is the area of the space, we are returning a %
  } 
})



svglite("figures/figure_S10.svg", width = 4.1, height = 5.5)

ggplot(data, aes(x=individuals, y=outcomeProportion, fill=outcomeGroup)) + 
      scale_x_continuous(breaks = round(seq(min(individuals), max(individuals), by = 15),1)) +
      xlab("Individuals") + 
      theme(panel.background = element_blank()) +
      theme(legend.position = "none") + 
      ylab("Outcome proportion") +
      scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
      geom_area() + 
      guides(fill=guide_legend(title="Outcome group"))+
      facet_grid(vars(fct_rev(loci),freq), vars(linked))
dev.off()