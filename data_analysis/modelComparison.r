
library(svglite)
library(ggplot2)
library(dplyr)


survivalW <- read.table("./data/final_23062025_model1.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise()
  
survivalSingle <- read.table("./data/final_23062025_model2.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise()

survivalRecurse <- read.table("./data/final_23062025_model3.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise()

slimData <- read.csv("./data/out_NumericalComparison_LIFPNG_1750435735.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise()

mc<-200
q<-0.04
gr<-1.5
i<-25
depFreq <- function(x) {filter(x, MutationCount==mc, Individuals == i, GrowthRate == gr)}
depLoci <- function(x) {filter(x, MutationFrequency==q, Individuals == i, GrowthRate == gr)}
depGrowth <- function(x) {filter(x, MutationCount==mc, MutationFrequency==q, Individuals == i)}
depInds <- function(x) {filter(x, MutationCount==mc, MutationFrequency==q, GrowthRate == gr)}

svglite("figures/figure_1.svg", width = 8, height = 6)
funcrange<-list(depFreq,c(0.0,0.1))
plot(extinctionRate~MutationFrequency,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab=expression("Deleterious recessive frequency (" * q * ")"),ylab="Proportion extinct")
dev.off()


svglite("figures/figure_S1A.svg", width = 8, height = 6)
funcrange<-list(depLoci,c(0,250))
plot(extinctionRate~MutationCount,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationCount,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationCount,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~MutationCount,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab=expression("Deleterious loci count (" * l * ")"),ylab="Proportion extinct")
dev.off()

svglite("figures/figure_S1B.svg", width = 8, height = 6)
funcrange<-list(depInds,c(0,120))
plot(extinctionRate~Individuals,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~Individuals,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~Individuals,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~Individuals,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab="Founding population size (N(0))",ylab="Proportion extinct")
dev.off()

svglite("figures/figure_S1C.svg", width = 8, height = 6)
funcrange<-list(depGrowth,c(1.0,2.0))
plot(extinctionRate~GrowthRate,data=funcrange[[1]](survivalW),col="red",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=funcrange[[1]](survivalSingle),col="green",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=funcrange[[1]](survivalRecurse),col="blue",ylim=c(0,1),xlim=funcrange[[2]],pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=funcrange[[1]](slimData),col="black",ylim=c(0,1),xlim=funcrange[[2]],xlab=expression("Reproductive output (" * b * ")"),ylab="Proportion extinct")
dev.off()


# Broader comparison of best Numerical vs Slim hermaphrodite

survivalRecurse <- read.table("./data/bestmodeldiff_02072025_2_model3.csv", header = FALSE, sep=",") %>%
  rename(MutationFrequency=V2, MutationCount=V1, Individuals=V3, GrowthRate=V4,Result=V5) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise() %>% 
  arrange(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count)
slimData <- read.csv("./data/out_NumericalComparisonWide_LIFPNG_merged_1750792216_1751383220_1751461342.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count) %>% 
  summarise() %>% 
  arrange(MutationFrequency, MutationCount, Individuals, GrowthRate, extinctionRate, count)


all(survivalRecurse$count==400)
all(slimData$count==400)
all(survivalRecurse$MutationFrequency == slimData$MutationFrequency)
length(survivalRecurse$extinctionRate) == length(slimData$extinctionRate)
slimData["modelDiffs"] = slimData$extinctionRate-survivalRecurse$extinctionRate

summary(slimData$modelDiffs)

filteredSlimData <- slimData %>% filter(extinctionRate+modelDiffs/2 > 0.025,extinctionRate+modelDiffs/2 < 0.975)

summary(filteredSlimData$modelDiffs)


dotpointscale = 0.75/1.35
svglite("figures/figure_S2.svg", width = 13*dotpointscale, height = 10*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
hist(filteredSlimData$modelDiffs,ylim=c(0,100),xlim=c(-0.5,0.2),breaks = seq(from=-0.5025, to=0.2025, by=0.005),main=NULL,xlab="Difference in extinction proportions between individual based and numerical model datapoints")
dev.off()
