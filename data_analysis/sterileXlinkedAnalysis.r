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
simulatedDataset2 <- "./data/out_parameterSpace_LIFP_1715053595.csv" 
simulatedDataset <- "./data/out_parameterSpace_LIFPNG_1748848797.csv "


paramaterSpace <- read.csv(simulatedDataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, extinctionRate, count) %>%
  filter(2/GrowthRate <=(1-MutationFrequency^2)^MutationCount) %>%
  arrange(extinctionRate) %>%
  summarise()


unfilteredXTestData <-  read.csv(simulatedDataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) 
  
  
xtestData <- unfilteredXTestData  %>%
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, extinctionRate, expectedExtinctionRate, count) %>% 
  summarise()

sum(xtestData$count)/nrow(unfilteredXTestData)

xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$extinctionRate-xtestData[xtestData$Xlinked==0,]$extinctionRate
xlinkedRisk <- xtestData[xtestData$Xlinked==1,]$extinctionRate/xtestData[xtestData$Xlinked==0,]$extinctionRate

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




svglite("figures/figure_S1B.svg", width = 10, height = 10)
hist(xlinkedRisk[xlinkedRisk<=1000],xlim=c(0.125,4),log='x',breaks = 10^seq(from=-0.73, to=0.55, by=0.01),ylim=c(0,6),main=NULL,xlab="Risk ratio of extinction given xlinked instead of autosomal",xaxt='n')
abline(v = 1,col='#f55b02')
ticks_at =c(2^seq(from=-2, to=8, by=1),400)
axis(side=1, at=ticks_at, labels=ticks_at)
dev.off()

svglite("figures/figure_3A-big-Risk.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(xlinkedRisk,ylim=c(0,200),xlim=c(-0.25,0.25),breaks = seq(from=-0.40125, to=0.40125, by=0.0025),main=NULL,xlab="Difference in extinction rates between xlinked and autosomal datapoints",axes=FALSE)
dev.off()

# filter(GrowthRate < 6) %>% #%>%

unfilteredSterileTestData <- read.csv(simulatedDataset) %>%  filter( (Sterile==0&PostCompetitionMutationTiming==1)|(Sterile==0&PostCompetitionMutationTiming==0))%>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) 
  
  
steriletestData <- unfilteredSterileTestData  %>%
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, expectedExtinctionRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, expectedExtinctionRate, extinctionRate, count) %>% 
  summarise()


sum(steriletestData$count)/nrow(unfilteredSterileTestData)

sterileDiffs <- steriletestData[steriletestData$PostCompetitionMutationTiming==1,]$extinctionRate-steriletestData[steriletestData$PostCompetitionMutationTiming==0,]$extinctionRate
# sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$extinctionRate-steriletestData[steriletestData$Sterile==0,]$extinctionRate
sterileRisk <- steriletestData[steriletestData$Sterile==1,]$extinctionRate/steriletestData[steriletestData$Sterile==0,]$extinctionRate
sterileLogitDiff <- log(steriletestData[steriletestData$Sterile==1,]$extinctionRate/(1-steriletestData[steriletestData$Sterile==1,]$extinctionRate)) - 
  log(steriletestData[steriletestData$Sterile==0,]$extinctionRate/(1-steriletestData[steriletestData$Sterile==0,]$extinctionRate))


nsamples <- 400
s<-(-nsamples):nsamples
p <- steriletestData[steriletestData$PostCompetitionMutationTiming==1,]$expectedExtinctionRate
y<-sapply(s, function(x) sum(diffBin(x, nsamples, p, nsamples, p)))
write.csv(data.frame(s,y),"./data/expectedNullPrePostLethalDiff.csv")

s <- read.csv('./data/expectedNullSterileLethalDiff.csv')$s
y <- read.csv('./data/expectedNullSterileLethalDiff.csv')$y
mean(sterileDiffs)
sterileLessData <- subset(steriletestData[steriletestData$Sterile==1,], select = -c(Sterile, expectedExtinctionRate))
sterileLessData <- subset(steriletestData[steriletestData$PostCompetitionMutationTiming==1,], select = -c(Sterile, expectedExtinctionRate))
cor(sterileDiffs,sterileLessData)
sum(y) / sum(table(sterileDiffs)) 

dotpointscale = 0.75/1.35
svglite("figures/figure_3B-olxel.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileDiffs,ylim=c(0,100),xlim=c(-0.1,1.0),breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xlab="Difference in extinction rates between postlethal and prelethal datapoints",axes=FALSE)
par(new=TRUE)
plot(s/nsamples,y,ylim=c(0,100),xlim=c(-0.1,1.0),xlab='',ylab='',xaxp=c(-0.2, 1.0, 12),cex=dotpointscale/4)
dev.off()

svglite("figures/figure_3B2-big.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileRisk[sterileRisk!=Inf],ylim=c(0,200),xlim=c(-0.5,100),breaks = seq(from=0, to=370, by=1),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints")
dev.off()
svglite("figures/figure_3B3-big.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileLogitDiff[sterileLogitDiff!=Inf],ylim=c(0,150),xlim=c(-1.0,9.0),breaks = seq(from=-1.00125, to=10.0125, by=0.05),main=NULL,xlab="Difference in logit extinction rates between sterile and lethal datapoints")
dev.off()

sterileRisk <- steriletestData[steriletestData$Sterile==1,]$extinctionRate/steriletestData[steriletestData$Sterile==0,]$extinctionRate
# hist(sterileRisk,ylim=c(0,600),xlim=c(-0.2,1.0),breaks = seq(from=-1.005, to=0.905, by=0.01),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)

cor(sterileRisk,sterileLessData)


svglite("figures/figure_S1.svg", width = 10, height = 10)
hist(sterileRisk[sterileRisk!=Inf],xlim=c(0.25,400),log='x',breaks = 10^seq(from=-0.4, to=2.56, by=0.01),ylim=c(0,2.5),main=NULL,xlab="Risk ratio of extinction given sterility instead of lethality",xaxt='n')
abline(v = 1,col='#f55b02')
ticks_at =c(2^seq(from=-2, to=8, by=1),400)
axis(side=1, at=ticks_at, labels=ticks_at)
dev.off()

hist(sterileRisk,breaks = seq(from=-1.005, to=0.905, by=0.05),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
10^mean(log10(sterileRisk[sterileRisk<=1000000]))
summary(sterileRisk[sterileRisk<=100000000])
