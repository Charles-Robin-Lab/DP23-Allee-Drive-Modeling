library(dplyr)

groupedData <- read.csv("./data/out.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()


# Extra stats
groupedData$logisticSurvivalRate <- log(groupedData$survivalRate/(1-groupedData$survivalRate))
groupedData$sQIndividuals <- groupedData$Individuals^2
groupedData$sQMutationFrequency <- groupedData$MutationFrequency^2
groupedData$eIndividuals <- exp(groupedData$Individuals)
groupedData$individualChance <- 0.5*(1-(1-groupedData$MutationFrequency^2)^groupedData$MutationCount)


#Model decission
# library(pROC)
# roc.curve <- roc()

options(width=2000)
model.slice <- glm(survivalRate~MutationCount+I(MutationCount^2),family=binomial(),data=dataSlice2d.mutationCoun)
model.all <- glm(survivalRate~MutationCount*MutationFrequency*Individuals*GrowthRate*Sterile,family=binomial(),data=groupedData)
model.all <- glm(survivalRate~individualChance*MutationCount*MutationFrequency*(I(Individuals^2)+Individuals)*GrowthRate*Sterile,family=binomial(),data=groupedData)
model.allaic <- step(model.all,direction='both') 

model.sig <- glm(survivalRate ~ 
  GrowthRate + 
  MutationCount:MutationFrequency + 
  MutationCount:GrowthRate +
  MutationFrequency:GrowthRate + 
  Individuals:GrowthRate + 
  GrowthRate:Sterile + 
  MutationCount:MutationFrequency:Individuals + 
  MutationCount:Individuals:GrowthRate  + 
  MutationFrequency:Individuals:GrowthRate  + 
  Individuals:GrowthRate:Sterile + 
  MutationCount:MutationFrequency:Individuals:GrowthRate + 
  MutationCount:MutationFrequency:Individuals:Sterile + 
  MutationCount:Individuals:GrowthRate:Sterile  + 
  MutationFrequency:Individuals:GrowthRate:Sterile )

model.all.added <-glm(survivalRate~MutationCount*(I(MutationFrequency^2)+MutationFrequency)*(I(Individuals^2)+Individuals)*GrowthRate*Sterile,family=binomial(),data=groupedData)


library(ggplot2)


# linked
dataset <-"./data/out_LoadEquivalenceParameterSpace_LIFP_1746635945.csv"
title <- "linked"
# unlinked
# dataset <-"./data/out_LoadEquivalenceParameterSpace_LIFP_1746978097.csv"
# title <- "unlinked"

groupedGraphData <- read.csv(dataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, extinctionRate, count) %>% 
  summarise()

groupedGraphData$extinctionRate <- ifelse(2/groupedGraphData$GrowthRate <= (1-groupedGraphData$MutationFrequency^2)^groupedGraphData$MutationCount,groupedGraphData$extinctionRate, -0.1)
groupedGraphData$extinctionRate <- ifelse(groupedGraphData$MutationFrequency==0.1 & groupedGraphData$MutationCount==400, 1.0, groupedGraphData$extinctionRate)

gr<-3
i<-25

focusedDataByLoad <- groupedGraphData %>%
filter(GrowthRate == gr) %>% 
filter(Individuals == i)
colors <- c("#ffffff", 
          "#ff0000", 
          "#00ff00",
          "#0000ff", 
          "#ffff00",  
          "#00ffff", 
          "#000000", 
          "#00ffff", 
          "#ffff00", 
          "#0000ff", 
          "#00ff00", 
          "#ff0000")

ggplot(focusedDataByLoad,aes(x=MutationFrequency,y=MutationCount),) +
ggtitle(sprintf("%s gr=%s i=%s",title,gr,i)) +
geom_tile(aes(fill=extinctionRate)) +
# scale_x_continuous(breaks = seq(0, max(focusedDataByLoad$Individuals), 25),expand = c(0,0)) +
# scale_y_continuous(breaks = seq(0, 150, 5) ,expand = c(0,0)) +
scale_fill_gradientn(name = "", colours = colors ,guide="none") +
theme_classic() +
labs(fill = "Generations") +
theme(legend.title.align=0.5)

print(arrange(middlesubset,extinctionRate)[c("extinctionRate","MutationCount","MutationFrequency")],n=1000)
diff <- 0.475

focusedDataByLoad$MeanIndividualMuts <- focusedDataByLoad$MutationFrequency *focusedDataByLoad$MutationCount
focusedDataByLoad$FilterRate <- 1-(1-focusedDataByLoad$MutationFrequency^2)^focusedDataByLoad$MutationCount
focusedDataByLoad$Hyperbola <- ((8.97+2/900)/focusedDataByLoad$MutationFrequency-40-80/9)/focusedDataByLoad$MutationCount
middlesubset <- filter(focusedDataByLoad,extinctionRate> 0.5-diff,extinctionRate< 0.5+diff)

length(unique(middlesubset$MutationCount))
length(middlesubset$MutationCount)
# length(unique(middlesubset$MutationFrequency))
model<-glm(extinctionRate~FilterRate*MeanIndividualMuts+MutationFrequency+MutationCount,family = binomial(),data=focusedDataByLoad)

step(model,direction='both')
model<-glm(extinctionRate~FilterRate+MeanIndividualMuts,family = binomial(),data=middlesubset)
# model<-glm(extinctionRate~FilterRate+MeanIndividualMuts+FilterRate:MeanIndividualMuts,family = binomial(),data=focusedDataByLoad)
summary(model)
