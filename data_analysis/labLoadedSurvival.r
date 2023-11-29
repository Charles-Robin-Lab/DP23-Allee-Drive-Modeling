if (!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales)

recombinationComparisonData <- read.csv("./data/out_LouisRecombinationRate_LLP_1699964632.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, CarryingCapacity, RecombinationRate, ChromosomeCount, MaxGenerations) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, CarryingCapacity, ChromosomeCount, MaxGenerations, survivalRate, loadSurvivalRate) %>% 
  summarise()


# Does the trans and cis data matter?
recombinationComparisonDataNonTrans <- read.csv("./data/out_LouisRecombinationRate_Nontrans_LLP_1700750219.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, CarryingCapacity, ChromosomeCount, MaxGenerations) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, CarryingCapacity, ChromosomeCount, MaxGenerations, survivalRate, loadSurvivalRate) %>% 
  summarise()



inds=2
inds=6
inds=10
transCompare <- subset(filter(recombinationComparisonDataNonTrans,Individuals==inds,MaxGenerations==2000), select = -c(MaxGenerations,Individuals) )
transCompare["survivalRate"] <- 100*transCompare["survivalRate"]
transCompare["loadSurvivalRate"] <- 100*transCompare["loadSurvivalRate"]
transCompare["transSurvivalIncrease"] <- 100*(filter(recombinationComparisonData,Individuals==inds,MaxGenerations==2000)["survivalRate"] - filter(recombinationComparisonDataNonTrans,Individuals==inds,MaxGenerations==2000)["survivalRate"])

transCompare["transLoadSurvivalIncrease"] <- 100*(filter(recombinationComparisonData,Individuals==inds,MaxGenerations==2000)["loadSurvivalRate"] - filter(recombinationComparisonDataNonTrans,Individuals==inds,MaxGenerations==2000)["loadSurvivalRate"])


print(transCompare,n=128)
print(transCompare[order(transCompare$transSurvivalIncrease),],n=64)
print(transCompare[order(transCompare$transLoadSurvivalIncrease),],n=64)
# Answer: Yes but it varries wildly. A key pattern is that linked datapoints show a large decrease in survival when switching to trans, corresponding unlinked datapoints dont see this as strongly

# Does the carrying Capacity have a big effect anywhere?
# Hypothesis: higher carrying capacity = less small population competition = slightly more survival 
CCsurvivaldiffs <- filter(recombinationComparisonData,CarryingCapacity==1000,MaxGenerations==2000)$survivalRate-filter(recombinationComparisonData,CarryingCapacity==100,MaxGenerations==2000)$survivalRate
par(mfrow = c(1,1))
hist(CCsurvivaldiffs,breaks = seq(from=-1.005, to=1.005, by=0.01))
print(filter(recombinationComparisonData,CarryingCapacity==1000,MaxGenerations==2000)[abs(CCsurvivaldiffs)>0.5,],n=100)
CCNonExtinctdiffs <- CCsurvivaldiffs + filter(recombinationComparisonData,CarryingCapacity==1000,MaxGenerations==2000)$loadSurvivalRate-filter(recombinationComparisonData,CarryingCapacity==100,MaxGenerations==2000)$loadSurvivalRate
hist(CCNonExtinctdiffs,breaks = seq(from=-1.005, to=1.005, by=0.01))
print(filter(recombinationComparisonData,CarryingCapacity==1000,MaxGenerations==2000)[abs(CCNonExtinctdiffs)>0.2,],n=100)
# Answer: Nearly everywhere this seems to be true, there are a few(~10) datapoints where extinction is largely different, mostly high individual counts and low growth rates, all linked. Makes sense I think because the carrying capacity


# Questions
# - how different is loaded survival in different maxgens
# - is the 
# - we do see a decrease in loaded survival when carrying capcity is lowered at low maxgens right?
# - what is the best growth rate (range?)(time in vial) and mutation spread(lower first) for different migrant sizes(start at 2)? Best being biggest difference in loaded and unloaded survival
# - do we see more extinction when spread and more loaded survival when grouped
filter(recombinationComparisonData,Individuals==2)
filter(recombinationComparisonData,Individuals==2,MaxGenerations==100)
# varry growthrate on x
# diff linked unlinked

# format individuals x loci


# - how different is loaded survival in different maxgens
# recombinationComparisonData['outcome'] = rep("EXTINCTION")
# recombinationComparisonData['outcome'] = rep("LOADED_SURVIVAL")
# recombinationComparisonData['outcome'] = rep("SURVIVAL")
# Stacked
ggplot(recombinationComparisonData, aes(fill=outcome, y=survivalRate, x=MaxGenerations)) + 
    geom_bar(position="stack", stat="identity")


par(mfrow = c(3, 5))
recombinationComparisonDataGen <- filter(recombinationComparisonData,MaxGenerations==2000)
for (individualCount in c(2,6,10)) {
   for (lociCount in 2:5) {
      
    lociData <- filter(recombinationComparisonDataGen,MutationCount==lociCount,Individuals==individualCount)
    linked <- filter(lociData,ChromosomeCount==1)
    unlinked <- filter(lociData,ChromosomeCount!=1)
    survivalRateDiffs <- linked$survivalRate - unlinked$survivalRate
    loadedSurvivalRateDiffs <- linked$loadSurvivalRate - unlinked$loadSurvivalRate
      extinctionRateDiffs <- (1 - linked$loadSurvivalRate - linked$survivalRate) - (1 - unlinked$loadSurvivalRate - unlinked$survivalRate)
    
  # create data
  growthRate <- rep(linked$GrowthRate,times=3)
  outcomeProportion <- c(extinctionRateDiffs,loadedSurvivalRateDiffs,survivalRateDiffs)
  outcomeGroup <- rep(c("EXTINCTION","LOADED_SURVIVAL","SURVIVAL"),each=length(linked$GrowthRate))
  data <- data.frame(growthRate, abs(outcomeProportion), outcomeGroup)


  # stacked area chart
  print(ggplot(data,log="x", aes(x=growthRate, y=abs(outcomeProportion), fill=outcomeGroup)) + 
      geom_area() +
      ggtitle(sprintf("Inds %d , loci %d",individualCount,lociCount)) +
      # scale_x_continuous(expand = c(0, 0)) +
      scale_x_continuous(trans='log2',limits = c(min(growthRate), max(growthRate)),
      breaks = trans_breaks("log2", function(x) 2^(x)),expand = c(0,0)) +
      scale_y_continuous(expand = c(0, 0)) +
      annotation_logticks(sides="b"))
  
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



