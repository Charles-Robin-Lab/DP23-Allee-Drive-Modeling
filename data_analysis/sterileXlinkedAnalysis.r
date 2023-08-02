library(dplyr)





xtestData <- read.csv("../data/out_49786075.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  filter(survivalRate<0.5 && survivalRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()

xtestData <- read.csv("../data/out_49786075.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  filter(survivalRate==0.5) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()


xtestData <- read.csv("./data/out_49786075.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(expectedSurvivalRate = sum(Result != "EXTINCT") / count) %>%
  filter(expectedSurvivalRate<0.975 && expectedSurvivalRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count,expectedSurvivalRate) %>% 
  summarise()

xlessData <- subset(xtestData[xtestData$Xlinked==1,], select = -c(Xlinked, survivalRate))
xlessData$xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$survivalRate-xtestData[xtestData$Xlinked==0,]$survivalRate
cor(xlessData$xlinkedDiffs,xlessData)
hist(xlessData$xlinkedDiffs,ylim=c(0,600),xlim=c(-0.40,0.40),breaks = seq(from=-0.405, to=0.405, by=0.01))

s<-(-100):100

p <- xtestData[xtestData$Xlinked==1,]$expectedSurvivalRate

y<-sapply(s, function(x) sum(diffBin(x, 100, p, 100, p)))
par(new=TRUE)
plot(s/100,y,ylim=c(0,600),xlim=c(-0.40,0.40))


steriletestData <- read.csv("./data/out_49786075.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(expectedSurvivalRate = sum(Result != "EXTINCT") / count) %>%
  filter(expectedSurvivalRate<0.975 && expectedSurvivalRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, expectedSurvivalRate) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result != "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate, expectedSurvivalRate, count) %>% 
  summarise()




sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$survivalRate-steriletestData[steriletestData$Sterile==0,]$survivalRate
hist(sterileDiffs,ylim=c(0,600),xlim=c(-1.0,0.40),breaks = seq(from=-1.005, to=0.905, by=0.01))
boxplot(sterileDiffs,xlessData$xlinkedDiffs,ylim=c(0,600),xlim=c(-0.80,0.40))




# https://gist.github.com/coppeliaMLA/9681819

diffBin<-function(z, n1, p1, n2, p2){

  prob<-0
  
  if (z>=0){  
    for (i in 0:n1){     
      prob<-prob+dbinom(i+z, n1, p1)*dbinom(i, n2, p2)
    }
    
  }
  else
  {
    for (i in 0:n2){     
      prob<-prob+dbinom(i+z, n1, p1)*dbinom(i, n2, p2)
    }
  }
  return(prob)
}

#Example
s<-(-100):100

p <- steriletestData[steriletestData$Sterile==1,]$expectedSurvivalRate

p<-sapply(s, function(x) sum(diffBin(x, 100, p, 100, p)))
par(new=TRUE)
plot(s/100,p,ylim=c(0,600),xlim=c(-1.0,0.40))


