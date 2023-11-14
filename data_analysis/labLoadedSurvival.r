if (!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(dplyr)

recombinationComparisonData <- read.csv("./data/out_LouisRecombinationRate_LLP_1699600504.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount, MaxGenerations) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, ChromosomeCount, MaxGenerations, survivalRate, loadSurvivalRate, count) %>% 
  summarise()





# Questions
# - how different is loaded survival in different maxgens
# - we do see a decrease in loaded survival when carrying capcity is lowered at low maxgens right?
# - what is the best growth rate (range?)(time in vial) and mutation spread(lower first) for different migrant sizes(start at 2)? Best being biggest difference in loaded and unloaded survival
# - do we see more extinction when spread and more loaded survival when grouped
filter(recombinationComparisonData,Individuals==2)
filter(recombinationComparisonData,Individuals==2,MaxGenerations==100)
# varry growthrate on x
# diff linked unlinked

# format individuals x loci

par(mfrow = c(3, 5))

for (individualCount in c(2,6,10)) {
   for (lociCount in 2:5) {
      
    lociData <- filter(recombinationComparisonData,MutationCount==lociCount,Individuals==individualCount)
    linked <- filter(lociData,ChromosomeCount==1)
    unlinked <- filter(lociData,ChromosomeCount!=1)
    survivalRateDiffs <- linked$survivalRate - unlinked$survivalRate
    loadedSurvivalRateDiffs <- linked$loadSurvivalRate - unlinked$loadSurvivalRate
      extinctionRateDiffs <- (1 - linked$loadSurvivalRate - linked$survivalRate) - (1 - unlinked$loadSurvivalRate - unlinked$survivalRate)
    
  # create data
  GrowthRate <- rep(linked$GrowthRate,times=3)
  outcomeProportion <- c(extinctionRateDiffs,loadedSurvivalRateDiffs,survivalRateDiffs)
  outcomeGroup <- rep(c("EXTINCTION","LOADED_SURVIVAL","SURVIVAL"),each=length(linked$GrowthRate))
  data <- data.frame(GrowthRate, outcomeProportion, outcomeGroup)


  # stacked area chart
  ggplot(data,log="x", aes(x=GrowthRate, y=outcomeProportion, fill=outcomeGroup)) + 
      geom_area() +
      scale_x_continuous(trans='log2',limits = c(min(GrowthRate), max(GrowthRate)),
      breaks = trans_breaks("log2", function(x) 2^(x)),expand = c(0,0)) +
      scale_y_continuous(expand = c(0, 0)) +
      annotation_logticks(sides="b")

   }
}


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

for (lociCount in 2:5) {
   lociData <- filter(recombinationComparisonData,MutationCount==lociCount)
   linked <- filter(lociData,ChromosomeCount==1)
   unlinked <- filter(lociData,ChromosomeCount!=1)


}

# # create a dataset
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)
 
# # Grouped
# ggplot(data, aes(fill=condition, y=value, x=specie)) + 
#     geom_bar(position="dodge", stat="identity")


# Packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales)
 
# create data
autosomalRecombinationComparisonData<- recombinationComparisonData
GrowthRate <- rep(autosomalRecombinationComparisonData$GrowthRate,times=3)
outcomeProportion <- c(1-autosomalRecombinationComparisonData$survivalRate - autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$survivalRate)
outcomeGroup <- rep(c("EXTINCTION","LOADED_SURVIVAL","SURVIVAL"),each=length(autosomalRecombinationComparisonData$GrowthRate))
data <- data.frame(GrowthRate, outcomeProportion, outcomeGroup)


# stacked area chart
ggplot(data,log="x", aes(x=GrowthRate, y=outcomeProportion, fill=outcomeGroup)) + 
    geom_area() +
    scale_x_continuous(trans='log2',limits = c(min(GrowthRate), max(GrowthRate)),
    breaks = trans_breaks("log2", function(x) 2^(x)),expand = c(0,0)) +
    scale_y_continuous(expand = c(0, 0)) +
    annotation_logticks(sides="b") +
    geom_vline(xintercept=2.1e-8, color="yellow")



