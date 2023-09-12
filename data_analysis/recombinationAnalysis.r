
library(dplyr)

recombinationComparisonData <- read.csv("./data/BHSRecomb100.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, survivalRate, loadSurvivalRate,count) %>% 
  summarise()




plot(survivalRate+loadSurvivalRate~RecombinationRate,ylab="Surivival Rate",data=filter(recombinationComparisonData, Xlinked==0),ylim=c(0,1),log="x")
par(new=TRUE)
plot(loadSurvivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 19,ylim=c(0,1),log="x")
par(new=TRUE)
plot(survivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 2,ylim=c(0,1),log="x")
abline(v=0.00101, col="blue")
abline(v=1.0e-5, col="red")

plot(survivalRate+loadedSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),ylim=c(0,1),log="x")
par(new=TRUE)
plot(loadSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 19,ylim=c(0,1),log="x")
par(new=TRUE)
plot(survivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 2,ylim=c(0,1),log="x")



# Packages
library(ggplot2)
library(dplyr)
library(scales)
 
# create data
autosomalRecombinationComparisonData<- filter(recombinationComparisonData, Xlinked==0)
RecombinationRate <- rep(autosomalRecombinationComparisonData$RecombinationRate,times=3)
outcomeProportion <- c(1-autosomalRecombinationComparisonData$survivalRate - autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$survivalRate)
outcomeGroup <- rep(c("EXTINCTION","LOADED_SURVIVAL","SURVIVAL"),each=length(autosomalRecombinationComparisonData$RecombinationRate))
data <- data.frame(RecombinationRate, viability, outcomeGroup)

# stacked area chart
ggplot(data,log="x", aes(x=RecombinationRate, y=outcomeProportion, fill=outcomeGroup)) + 
    geom_area() +
    scale_x_continuous(trans='log10',limits = c(min(RecombinationRate[RecombinationRate!=0]), max(RecombinationRate)),
    breaks = trans_breaks("log10", function(x) 10^(x)),expand = c(0,0)) +
    scale_y_continuous(expand = c(0, 0)) +
    annotation_logticks(sides="b")
