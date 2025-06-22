
library(svglite)
library(ggplot2)
library(dplyr)


survivalW <- read.table("./data/finalish1.5_21062025_model1.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionProbability = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionProbability, count) %>% 
  summarise()
  
survivalSingle <- read.table("./data/finalish1.5_21062025_model2.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionProbability = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionProbability, count) %>% 
  summarise()


survivalRecurse <- read.table("./data/finalish1.5_21062025_model3.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionProbability = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionProbability, count) %>% 
  summarise()

slimData <- read.csv("./data/out_AnalyticalComparison_LIFPNG_1750435735.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionProbability = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionProbability, count) %>% 
  summarise()

mc<-200
q<-0.04
gr<-1.5
i<-25
depFreq <- function(x) {filter(x, MutationCount==mc, Individuals == i, GrowthRate == gr)}
depLoci <- function(x) {filter(x, MutationFrequency==q, Individuals == i, GrowthRate == gr)}
depGrowth <- function(x) {filter(x, MutationCount==mc, MutationFrequency==q, Individuals == i)}
depInds <- function(x) {filter(x, MutationCount==mc, MutationFrequency==q, GrowthRate == gr)}

# dataSlice2d.MutationFrequency <- read.csv("./data/modelslim.csv")
# A = data.frame(x = freqRange, y = survivalW)
# B = data.frame(x = freqRange, y = survivalSingle)
# C = data.frame(x = freqRange, y = survivalRecurse)
# D = data.frame(x = dataSlice2d.MutationFrequency$MutationFrequency, y = dataSlice2d.MutationFrequency$extinctionProbability)

# A = data.frame(x = rnorm(10),y=rnorm(10))
# B = data.frame(x = rnorm(10),y=rnorm(10))
# ggplot(A,aes(x,y)) +geom_point() +geom_point(data=B,colour='red') + xlim(0, 10)







# freqRange <- slimData$MutationFrequency
# svglite("figures/figure_1.2_1747151880.svg", width = 8, height = 6)
# dev.off()

funcrange<-list(depFreq,c(0.0,0.1))

plot(extinctionProbability~MutationFrequency,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationFrequency,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationFrequency,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationFrequency,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab="Deleterious recessive frequency",ylab="Extinction probability")


funcrange<-list(depLoci,c(0,250))
plot(extinctionProbability~MutationCount,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationCount,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationCount,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationCount,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab="Loci Count",ylab="Extinction probability")


funcrange<-list(depInds,c(0,120))
plot(extinctionProbability~Individuals,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~Individuals,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~Individuals,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~Individuals,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab="Individual count",ylab="Extinction probability")


funcrange<-list(depGrowth,c(1.0,2.0))
plot(extinctionProbability~GrowthRate,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~GrowthRate,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~GrowthRate,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~GrowthRate,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab="GrowthRate",ylab="Extinction probability")
