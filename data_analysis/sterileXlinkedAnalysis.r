library(dplyr)
library(car)
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



xtestData <- read.csv("./data/out_50022325_2.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) %>%
  filter(expectedExtinctionRate<0.975 && expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, extinctionRate,expectedExtinctionRate, count) %>% 
  summarise()

xlessData <- subset(xtestData[xtestData$Xlinked==1,], select = -c(Xlinked, expectedExtinctionRate))
xlessData$xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$extinctionRate-xtestData[xtestData$Xlinked==0,]$extinctionRate
cor(xlessData$xlinkedDiffs,xlessData)
par(mar=c(5,4,2,2)+0.1)
hist(xlessData$xlinkedDiffs,ylim=c(0,600),xlim=c(-0.25,0.25),breaks = seq(from=-0.405, to=0.405, by=0.01),main=NULL,xlab="Difference in extinction rates between xlinked and autosomal datapoints",axes=FALSE)

s1<-(-100):100

p1 <- xtestData[xtestData$Xlinked==1,]$expectedExtinctionRate
# n1 <- xtestData[xtestData$Xlinked==1,]$count
# n2 <- xtestData[xtestData$Xlinked==0,]$count
# y  <-rep(0, length(p))
# for (i in seq_along(p)) {
#   y[i] <-sapply(s, function(x) sum(diffBin(x, n1[i], p[i], n2[i], p[i])))
# }
y1<-sapply(s1, function(x) sum(diffBin(x, 100, p1, 100, p1)))
par(new=TRUE)
plot(s1/100,y1,ylim=c(0,600),xlim=c(-0.25,0.25),xlab='',ylab='',xaxp=c(-0.25, 0.25, 10),cex=1.0)


steriletestData <- read.csv("./data/out_50022325_2.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) %>%
  filter(expectedExtinctionRate<0.975 && expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, expectedExtinctionRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, expectedExtinctionRate, extinctionRate, count) %>% 
  summarise()




sterileLessData <- subset(steriletestData[steriletestData$Sterile==1,], select = -c(Sterile, expectedExtinctionRate))
sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$extinctionRate-steriletestData[steriletestData$Sterile==0,]$extinctionRate
cor(sterileDiffs,sterileLessData)
hist(sterileDiffs,ylim=c(0,600),xlim=c(-0.2,1.0),breaks = seq(from=-1.005, to=0.905, by=0.01),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
# boxplot(sterileDiffs,xlessData$xlinkedDiffs,ylim=c(0,600),xlim=c(-0.80,0.40))

sterileLessData <- subset(steriletestData[steriletestData$Sterile==1,], select = -c(Sterile, expectedExtinctionRate))
sterileRisk <- steriletestData[steriletestData$Sterile==1,]$extinctionRate/steriletestData[steriletestData$Sterile==0,]$extinctionRate
cor(sterileRisk,sterileLessData)
hist(sterileRisk,ylim=c(0,600),xlim=c(-0.2,1.0),breaks = seq(from=-1.005, to=0.905, by=0.01),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)


hist(sterileRisk[sterileRisk<=1000],xlim=c(0.25,400),log='x',breaks = 10^seq(from=-0.4, to=2.5, by=0.01),main=NULL,xlab="Risk ratio of extinction given sterility instead of lethality",xaxt='n')
abline(v = 1,col='#f55b02')
axis(side=1, at=2^seq(from=-2, to=8, by=1), labels=2^seq(from=-2, to=8, by=1))

hist(sterileRisk,breaks = seq(from=-1.005, to=0.905, by=0.05),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
10^mean(log10(sterileRisk[sterileRisk<=1000000]))
summary(sterileRisk[sterileRisk<=100000000])

#Example
s<-(-100):100

p <- steriletestData[steriletestData$Sterile==1,]$expectedExtinctionRate

p<-sapply(s, function(x) sum(diffBin(x, 100, p, 100, p)))
par(new=TRUE)
plot(s/100,p,ylim=c(0,600),xlim=c(-0.2,1.0),xlab='',ylab='',xaxp=c(-0.2, 1.0, 12),cex=0.5)


s<-(-400):400

p <- steriletestData[steriletestData$Sterile==1,]$expectedExtinctionRate

p<-sapply(s, function(x) sum(diffBin(x, 400, p, 400, p)))
par(new=TRUE)
plot(s/400,p,ylim=c(0,600),xlim=c(-0.2,1.0),xlab='',ylab='',xaxp=c(-0.2, 1.0, 12),cex=0.5)

