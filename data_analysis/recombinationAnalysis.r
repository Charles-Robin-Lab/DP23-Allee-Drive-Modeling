
recombinationComparisonData <- read.csv("./data/recombinationComparison2.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED" ) / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, survivalRate, loadSurvivalRate,count) %>% 
  summarise()




plot(survivalRate+loadedSurvivalRate~RecombinationRate,ylab="Surivival Rate",data=filter(recombinationComparisonData, Xlinked==0),ylim=c(0,1),log="x")
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
viability <- c(1-autosomalRecombinationComparisonData$survivalRate - autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$loadSurvivalRate,autosomalRecombinationComparisonData$survivalRate)
outcomeGroup <- rep(c("EXTINCTION","LOADED_SURVIVAL","SURVIVAL"),each=length(autosomalRecombinationComparisonData$RecombinationRate))
data <- data.frame(RecombinationRate, viability, outcomeGroup)

# stacked area chart
ggplot(data,log="x", aes(x=RecombinationRate, y=viability, fill=outcomeGroup)) + 
    geom_area() +
    scale_x_continuous(trans='log10',limits = c(1e-15, 1e-3),
    breaks = trans_breaks("log10", function(x) 10^(x)))
    + annotation_logticks()
