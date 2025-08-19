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
simulatedDataset <- "./data/out_bestParameterSpace_LIFPNG_1752545194.csv"
dataset <- read.csv(simulatedDataset)


unfilteredXTestData <-  dataset %>% 
  filter(PostCompetitionMutationTiming==0) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, PostCompetitionMutationTiming) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) 
  
  
xtestData <- unfilteredXTestData  %>%
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, extinctionRate, expectedExtinctionRate, count) %>% 
  summarise()
# % of unfiltered data
sum(xtestData$count)/nrow(unfilteredXTestData)
# roughly verify the order
all(xtestData[xtestData$Xlinked==1,]$MutationFrequency == xtestData[xtestData$Xlinked==0,]$MutationFrequency)

xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$extinctionRate-xtestData[xtestData$Xlinked==0,]$extinctionRate
xlinkedRisk <- xtestData[xtestData$Xlinked==1,]$extinctionRate/xtestData[xtestData$Xlinked==0,]$extinctionRate

xislowertail <- xtestData[xtestData$Xlinked==1,]$extinctionRate <= xtestData[xtestData$Xlinked==1,]$expectedExtinctionRate

nsamples <-400

a <- xtestData[xtestData$Xlinked==1,]$extinctionRate*nsamples
b <- xtestData[xtestData$Xlinked==0,]$extinctionRate*nsamples
c <- nsamples-a
d <- nsamples-b

p <- mapply(function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher.test(mat)$p.value
}, a, b, c, d)

pthresh <- 0.05/length(p)
pthresh <- 0.05
pthresh <- 1-(1-0.05)^(1/length(p))
print(pthresh)
xlower <- xtestData[xtestData$Xlinked==1,][p <pthresh & xislowertail,]
xhigher <- xtestData[xtestData$Xlinked==1,][p <pthresh & !xislowertail,]
print(colMeans(xtestData), digits = 4)
print(colMeans(xlower), digits = 4)
print(colMeans(xhigher), digits = 4)


if (!file.exists("./data/expectedNullXlinkedAutosomalDiff.csv")) {
  nsamples <- 400
  s1<-(-nsamples):nsamples
  p1 <- xtestData[xtestData$Xlinked==1,]$expectedExtinctionRate
  y1<-sapply(s1, function(x) sum(diffBin(x, nsamples, p1, nsamples, p1)))
  write.csv(data.frame(s1,y1),"./data/expectedNullXlinkedAutosomalDiff.csv")
}

s1 <- read.csv('./data/expectedNullXlinkedAutosomalDiff.csv')$s1
y1 <- read.csv('./data/expectedNullXlinkedAutosomalDiff.csv')$y1

pthresh <- 1-(1-0.05)^(1/(2*nsamples))
print(pthresh)
qcoltop<- mapply(function(p) {
  qbinom(1-pthresh, length(xlinkedDiffs), p) 
}, y1/length(xlinkedDiffs))
qcolbot<- mapply(function(p) {
  qbinom(pthresh, length(xlinkedDiffs), p) 
}, y1/length(xlinkedDiffs))

xlessData <- subset(xtestData[xtestData$Xlinked==1,], select = -c(Xlinked, expectedExtinctionRate))
cor(xlinkedDiffs,xlessData)
# check same size(within a float error tolerence)
sum(y1) - as.double(length(xlinkedDiffs)) < 1e-10
dotpointscale = 0.75/1.35
svglite("figures/figure_5A.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(xlinkedDiffs,ylim=c(0,1000),xlim=c(-0.275,0.275),breaks = seq(from=-0.40125, to=0.40125, by=0.0025),main=NULL,xlab="Difference in extinction rates between xlinked and autosomal datapoints",axes=FALSE)
par(new=TRUE)
arrows((s1/nsamples)[qcoltop!=0], qcolbot[qcoltop!=0], (s1/nsamples)[qcoltop!=0], qcoltop[qcoltop!=0],col="#665757", lwd =0.75, length=0.010, angle=90, code=3)
par(new=TRUE)
hist(2*(xlower$extinctionRate - xlower$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.275,0.275),col="red",breaks = seq(from=-0.40125, to=0.40125, by=0.0025),main=NULL,axes=FALSE,xlab='',ylab='')
par(new=TRUE)
hist(2*(xhigher$extinctionRate - xhigher$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.275,0.275),col="red",breaks = seq(from=-0.40125, to=0.40125, by=0.0025),main=NULL,xaxt="n",xlab='',ylab='')
ticks <- seq(-0.27, 0.27, by = 0.01)
axis(1, at = ticks, labels = FALSE, tck = -0.01) 
labels <- seq(-0.25, 0.25, by = 0.05)
axis(1, at = labels, labels = labels, tck = -0.02)
dev.off()



unfilteredSterileTestData <- dataset %>%  
  filter(PostCompetitionMutationTiming==0) %>%
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

sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$extinctionRate-steriletestData[steriletestData$Sterile==0,]$extinctionRate
sterileRisk <- steriletestData[steriletestData$Sterile==1,]$extinctionRate/steriletestData[steriletestData$Sterile==0,]$extinctionRate
sterileLogitDiff <- log(steriletestData[steriletestData$Sterile==1,]$extinctionRate/(1-steriletestData[steriletestData$Sterile==1,]$extinctionRate)) - 
  log(steriletestData[steriletestData$Sterile==0,]$extinctionRate/(1-steriletestData[steriletestData$Sterile==0,]$extinctionRate))


sterileislowertail <- steriletestData[steriletestData$Sterile==1,]$extinctionRate <= steriletestData[steriletestData$Sterile==1,]$expectedExtinctionRate


a <- steriletestData[steriletestData$Sterile==1,]$extinctionRate*nsamples
b <- steriletestData[steriletestData$Sterile==0,]$extinctionRate*nsamples
c <- nsamples-a
d <- nsamples-b

p <- mapply(function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher.test(mat)$p.value
}, a, b, c, d)

pthresh <- 0.05
pthresh <- 0.05/length(p)
pthresh <- 1-(1-0.05)^(1/length(p))
print(pthresh)
sterilelower <- steriletestData[steriletestData$Sterile==1,][p <pthresh & sterileislowertail,]
sterilehigher <- steriletestData[steriletestData$Sterile==1,][p <pthresh & !sterileislowertail,]


if (!file.exists("./data/expectedNullSterileLethalDiff.csv")) {
  nsamples <- 400
  s<-(-nsamples):nsamples
  p <- steriletestData[steriletestData$Sterile==1,]$expectedExtinctionRate
  y<-sapply(s, function(x) sum(diffBin(x, nsamples, p, nsamples, p)))
  write.csv(data.frame(s,y),"./data/expectedNullSterileLethalDiff.csv")
}

s <- read.csv('./data/expectedNullSterileLethalDiff.csv')$s
y <- read.csv('./data/expectedNullSterileLethalDiff.csv')$y
mean(sterileDiffs)
sterileLessData <- subset(steriletestData[steriletestData$Sterile==1,], select = -c(Sterile, expectedExtinctionRate))
cor(sterileDiffs,sterileLessData)
sum(y) / sum(table(sterileDiffs)) 


pthresh <- 1-(1-0.05)^(1/(2*nsamples))
print(pthresh)
qcoltop<- mapply(function(p) {
  qbinom(1-pthresh, length(sterileDiffs), p) 
}, y/length(sterileDiffs))
qcolbot<- mapply(function(p) {
  qbinom(pthresh, length(sterileDiffs), p) 
}, y/length(sterileDiffs))


dotpointscale = 0.75/1.35
svglite("figures/figure_5B.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileDiffs,ylim=c(0,1000),xlim=c(-0.1,1.0),breaks = seq(from=-1.00125, to=1.00125, by=0.0025),lwd =0.1,main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
par(new=TRUE)
hist(2*(sterilelower$extinctionRate - sterilelower$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,1.0),col="red",lwd =0.1,breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,axes=FALSE,xlab='',ylab='')
par(new=TRUE)
hist(2*(sterilehigher$extinctionRate - sterilehigher$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,1.0),col="red",breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xaxp=c(-0.2, 1.0, 12),xlab='',ylab='')
arrows((s/nsamples)[qcoltop!=0], qcolbot[qcoltop!=0], (s/nsamples)[qcoltop!=0], qcoltop[qcoltop!=0],col="#665757", lwd =0.75, length=0.010, angle=90, code=3)
dev.off()

cor(sterileRisk,sterileLessData)
range(sterileRisk[sterileRisk!=Inf])

svglite("figures/figure_S4.svg", width = 10, height = 10)
hist(sterileRisk[sterileRisk!=Inf],xlim=c(0.25,400),log='x',breaks = 10^seq(from=-0.6, to=2.6, by=0.01),ylim=c(0,2.5),main=NULL,xlab="Risk ratio of extinction given sterility instead of lethality",xaxt='n')
abline(v = 1,col='#f55b02')
ticks_at =c(2^seq(from=-2, to=8, by=1),400)
axis(side=1, at=ticks_at, labels=ticks_at)
dev.off()

hist(sterileRisk,breaks = seq(from=-1.005, to=0.905, by=0.05),main=NULL,xlab="Difference in extinction rates between sterile and lethal datapoints",axes=FALSE)
10^mean(log10(sterileRisk[sterileRisk<=1000000]))
summary(sterileRisk[sterileRisk<=100000000])

# pre and post diffs

unfilteredTimingTestData <- dataset %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(expectedExtinctionRate = sum(Result == "EXTINCT") / count) 
  
  
timingTestData <- unfilteredTimingTestData  %>% 
  filter(expectedExtinctionRate<0.975,expectedExtinctionRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, expectedExtinctionRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, PostCompetitionMutationTiming, expectedExtinctionRate, extinctionRate, count) %>% 
  summarise()


sum(timingTestData$count)/nrow(unfilteredTimingTestData)
1-sum(timingTestData$count[timingTestData$Sterile==TRUE])/nrow(unfilteredTimingTestData[unfilteredTimingTestData$Sterile==TRUE,])
1-sum(timingTestData$count[timingTestData$Sterile==FALSE])/nrow(unfilteredTimingTestData[unfilteredTimingTestData$Sterile==FALSE,])

lethalTimingDiffs <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,]$extinctionRate-timingTestData[timingTestData$PostCompetitionMutationTiming==0  & timingTestData$Sterile==0,]$extinctionRate

lethalTimingislowertail <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,]$extinctionRate <= timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,]$expectedExtinctionRate


a <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,]$extinctionRate*nsamples
b <- timingTestData[timingTestData$PostCompetitionMutationTiming==0 & timingTestData$Sterile==0,]$extinctionRate*nsamples
c <- nsamples-a
d <- nsamples-b

p <- mapply(function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher.test(mat)$p.value
}, a, b, c, d)

pthresh <- 0.05
pthresh <- 0.05/length(p)
pthresh <- 1-(1-0.05)^(1/length(p))
print(pthresh)
lethalTiminglower <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,][p <pthresh & lethalTimingislowertail,]
lethalTiminghigher <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0,][p <pthresh & !lethalTimingislowertail,]


sterileTimingDiffs <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,]$extinctionRate-timingTestData[timingTestData$PostCompetitionMutationTiming==0 & timingTestData$Sterile==1,]$extinctionRate

sterileTimingislowertail <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,]$extinctionRate <= timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,]$expectedExtinctionRate


a <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,]$extinctionRate*nsamples
b <- timingTestData[timingTestData$PostCompetitionMutationTiming==0 & timingTestData$Sterile==1,]$extinctionRate*nsamples
c <- nsamples-a
d <- nsamples-b

p <- mapply(function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher.test(mat)$p.value
}, a, b, c, d)

pthresh <- 0.05
pthresh <- 0.05/length(p)
pthresh <- 1-(1-0.05)^(1/length(p))
print(pthresh)
sterileTiminglower <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,][p <pthresh & sterileTimingislowertail,]
sterileTiminghigher <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,][p <pthresh & !sterileTimingislowertail,]


all((timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1,]$GrowthRate-timingTestData[timingTestData$PostCompetitionMutationTiming==0 & timingTestData$Sterile==1,]$GrowthRate) ==0)


if (!file.exists("./data/expectedNullPrePostLethalDiff.csv")) {
  nsamples <- 400
  s<-(-nsamples):nsamples
  p <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==0 ,]$expectedExtinctionRate
  y<-sapply(s, function(x) sum(diffBin(x, nsamples, p, nsamples, p)))
  write.csv(data.frame(s,y),"./data/expectedNullPrePostLethalDiff.csv")
}

s <- read.csv('./data/expectedNullPrePostLethalDiff.csv')$s
y <- read.csv('./data/expectedNullPrePostLethalDiff.csv')$y


pthresh <- 1-(1-0.05)^(1/(2*nsamples))
print(pthresh)
qcoltop<- mapply(function(p) {
  qbinom(1-pthresh, length(lethalTimingDiffs), p) 
}, y/length(lethalTimingDiffs))
qcolbot<- mapply(function(p) {
  qbinom(pthresh, length(lethalTimingDiffs), p) 
}, y/length(lethalTimingDiffs))

dotpointscale = 0.75/1.35
svglite("figures/figure_S9A.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(lethalTimingDiffs,ylim=c(0,1000),xlim=c(-0.1,0.4),breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xlab="Difference in extinction rates between post and premortality lethal datapoints",axes=FALSE)
par(new=TRUE)
hist(2*(lethalTiminglower$extinctionRate - lethalTiminglower$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,0.4),col="red",lwd =0.1,breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,axes=FALSE,xlab='',ylab='')
par(new=TRUE)
hist(2*(lethalTiminghigher$extinctionRate - lethalTiminghigher$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,0.4),col="red",breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xaxp=c(-0.2, 1.0, 12),xlab='',ylab='')
arrows((s/nsamples)[qcoltop!=0], qcolbot[qcoltop!=0], (s/nsamples)[qcoltop!=0], qcoltop[qcoltop!=0],col="#665757", lwd =0.75, length=0.010, angle=90, code=3)
dev.off()


if (!file.exists("./data/expectedNullPrePostSterileDiff.csv")) {
  nsamples <- 400
  s<-(-nsamples):nsamples
  p <- timingTestData[timingTestData$PostCompetitionMutationTiming==1 & timingTestData$Sterile==1 ,]$expectedExtinctionRate
  y<-sapply(s, function(x) sum(diffBin(x, nsamples, p, nsamples, p)))
  write.csv(data.frame(s,y),"./data/expectedNullPrePostSterileDiff.csv")
}


s <- read.csv('./data/expectedNullPrePostSterileDiff.csv')$s
y <- read.csv('./data/expectedNullPrePostSterileDiff.csv')$y

pthresh <- 1-(1-0.05)^(1/(2*nsamples))
print(pthresh)
qcoltop<- mapply(function(p) {
  qbinom(1-pthresh, length(sterileTimingDiffs), p) 
}, y/length(sterileTimingDiffs))
qcolbot<- mapply(function(p) {
  qbinom(pthresh, length(sterileTimingDiffs), p) 
}, y/length(sterileTimingDiffs))

dotpointscale = 0.75/1.35
svglite("figures/figure_S9B.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(sterileTimingDiffs,ylim=c(0,1000),xlim=c(-0.1,1.0),breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xlab="Difference in extinction rates between post and premortality sterility datapoints",axes=FALSE)
par(new=TRUE)
hist(2*(sterileTiminglower$extinctionRate - sterileTiminglower$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,1.0),col="red",lwd =0.1,breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,axes=FALSE,xlab='',ylab='')
par(new=TRUE)
hist(2*(sterileTiminghigher$extinctionRate - sterileTiminghigher$expectedExtinctionRate),ylim=c(0,1000),xlim=c(-0.1,1.0),col="red",breaks = seq(from=-1.00125, to=1.00125, by=0.0025),main=NULL,xaxp=c(-0.2, 1.0, 12),xlab='',ylab='')
arrows((s/nsamples)[qcoltop!=0], qcolbot[qcoltop!=0], (s/nsamples)[qcoltop!=0], qcoltop[qcoltop!=0],col="#665757", lwd =0.75, length=0.010, angle=90, code=3)
dev.off()
