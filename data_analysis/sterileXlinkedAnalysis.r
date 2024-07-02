library(svglite)
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

simulatedDataset <- "./data/out_parameterSpace_LIFP_1715053595.csv"

xtestData <- read.csv(simulatedDataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) %>%
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, extinctionRate,expectedExtinctionRate, count) %>% 
  summarise()

xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$extinctionRate-xtestData[xtestData$Xlinked==0,]$extinctionRate

nsamples <- 400
# s1<-(-nsamples):nsamples
# p1 <- xtestData[xtestData$Xlinked==1,]$expectedExtinctionRate
# y1<-sapply(s1, function(x) sum(diffBin(x, nsamples, p1, nsamples, p1)))
# write.csv(data.frame(s1,y1),"./data/expectedNullXlinkedAutosomalDiff.csv")

s1 <- read.csv('./data/expectedNullXlinkedAutosomalDiff.csv')$s1
y1 <- read.csv('./data/expectedNullXlinkedAutosomalDiff.csv')$y1

sum(y1) / sum(table(xlinkedDiffs))
xlessData <- subset(xtestData[xtestData$Xlinked==1,], select = -c(Xlinked, expectedExtinctionRate))
cor(xlinkedDiffs,xlessData)

dotpointscale = 0.75/1.35
svglite("figures/figure_3A-big.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(xlinkedDiffs,ylim=c(0,200),xlim=c(-0.25,0.25),breaks = seq(from=-0.40125, to=0.40125, by=0.0025),main=NULL,xlab="Difference in extinction rates between xlinked and autosomal datapoints",axes=FALSE)
par(new=TRUE)
plot(s1/nsamples,y1,ylim=c(0,200),xlim=c(-0.25,0.25),xlab='',ylab='',xaxp=c(-0.25, 0.25, 10),cex=dotpointscale/2)
dev.off()

steriletestData <- read.csv(simulatedDataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) %>%
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, expectedExtinctionRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, expectedExtinctionRate, extinctionRate, count) %>% 
  summarise()
sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$extinctionRate-steriletestData[steriletestData$Sterile==0,]$extinctionRate


# nsamples <- 400
# s<-(-nsamples):nsamples
# p <- steriletestData[steriletestData$Sterile==1,]$expectedExtinctionRate
# y<-sapply(s, function(x) sum(diffBin(x, nsamples, p, nsamples, p)))
# write.csv(data.frame(s,y),"./data/expectedNullSterileLethalDiff.csv")

s <- read.csv('./data/expectedNullSterileLethalDiff.csv')$s
y <- read.csv('./data/expectedNullSterileLethalDiff.csv')$y

sterileLessData <- subset(steriletestData[steriletestData$Sterile==1,], select = -c(Sterile, expectedExtinctionRate))
cor(sterileDiffs,sterileLessData)
sum(y) / sum(table(sterileDiffs)) 

dotpointscale = 0.75/1.35
svglite("figures/figure_3B-big.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileDiffs,ylim=c(0,200),xlim=c(-0.1,1.0),breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
par(new=TRUE)
plot(s/nsamples,y,ylim=c(0,200),xlim=c(-0.1,1.0),xlab='',ylab='',xaxp=c(-0.2, 1.0, 12),cex=dotpointscale/4)
dev.off()

sterileRisk <- steriletestData[steriletestData$Sterile==1,]$extinctionRate/steriletestData[steriletestData$Sterile==0,]$extinctionRate
# hist(sterileRisk,ylim=c(0,600),xlim=c(-0.2,1.0),breaks = seq(from=-1.005, to=0.905, by=0.01),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)

cor(sterileRisk,sterileLessData)


svglite("figures/figure_S1.svg", width = 10, height = 10)
hist(sterileRisk[sterileRisk<=1000],xlim=c(0.25,400),log='x',breaks = 10^seq(from=-0.4, to=2.56, by=0.01),ylim=c(0,2.5),main=NULL,xlab="Risk ratio of extinction given sterility instead of lethality",xaxt='n')
abline(v = 1,col='#f55b02')
ticks_at =c(2^seq(from=-2, to=8, by=1),400)
axis(side=1, at=ticks_at, labels=ticks_at)
dev.off()

hist(sterileRisk,breaks = seq(from=-1.005, to=0.905, by=0.05),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
10^mean(log10(sterileRisk[sterileRisk<=1000000]))
summary(sterileRisk[sterileRisk<=100000000])
