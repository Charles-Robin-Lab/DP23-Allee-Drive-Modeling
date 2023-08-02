
recombinationComparisonData <- read.csv("../data/recombinationComparison.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(anySurvivalRate = sum(Result == "SURVIVED" | Result == "LOADED_SURVIVAL") / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, anySurvivalRate, loadSurvivalRate,count) %>% 
  summarise()




plot(anySurvivalRate~RecombinationRate,ylab="Surivival Rate",data=filter(recombinationComparisonData, Xlinked==0),ylim=c(0,1),log="x")
par(new=TRUE)
plot(loadSurvivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 19,ylim=c(0,1),log="x")
par(new=TRUE)
plot(anySurvivalRate-loadSurvivalRate~RecombinationRate,ylab="",data=filter(recombinationComparisonData, Xlinked==0),pch = 2,ylim=c(0,1),log="x")
abline(v=0.00101, col="blue")
abline(v=1.0e-5, col="red")




plot(anySurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),ylim=c(0,1),log="x")
par(new=TRUE)
plot(loadSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 19,ylim=c(0,1),log="x")
par(new=TRUE)
plot(anySurvivalRate-loadSurvivalRate~RecombinationRate,data=filter(recombinationComparisonData, Xlinked==1),pch = 2,ylim=c(0,1),log="x")
  