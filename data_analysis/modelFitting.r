library(dplyr)
library(svglite)
library("gridExtra")

# groupedData <- read.csv("./data/out.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
#   mutate(count = n()) %>%
#   mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
#   summarise()


# # Extra stats
# groupedData$logisticSurvivalRate <- log(groupedData$survivalRate/(1-groupedData$survivalRate))
# groupedData$sQIndividuals <- groupedData$Individuals^2
# groupedData$sQMutationFrequency <- groupedData$MutationFrequency^2
# groupedData$eIndividuals <- exp(groupedData$Individuals)
# groupedData$individualChance <- 0.5*(1-(1-groupedData$MutationFrequency^2)^groupedData$MutationCount)


#Model decission
# library(pROC)
# roc.curve <- roc()

options(width=2000)
# model.slice <- glm(survivalRate~MutationCount+I(MutationCount^2),family=binomial(),data=dataSlice2d.mutationCoun)
# model.all <- glm(survivalRate~MutationCount*MutationFrequency*Individuals*GrowthRate*Sterile,family=binomial(),data=groupedData)
# model.all <- glm(survivalRate~individualChance*MutationCount*MutationFrequency*(I(Individuals^2)+Individuals)*GrowthRate*Sterile,family=binomial(),data=groupedData)
# model.allaic <- step(model.all,direction='both') 

# model.sig <- glm(survivalRate ~ 
#   GrowthRate + 
#   MutationCount:MutationFrequency + 
#   MutationCount:GrowthRate +
#   MutationFrequency:GrowthRate + 
#   Individuals:GrowthRate + 
#   GrowthRate:Sterile + 
#   MutationCount:MutationFrequency:Individuals + 
#   MutationCount:Individuals:GrowthRate  + 
#   MutationFrequency:Individuals:GrowthRate  + 
#   Individuals:GrowthRate:Sterile + 
#   MutationCount:MutationFrequency:Individuals:GrowthRate + 
#   MutationCount:MutationFrequency:Individuals:Sterile + 
#   MutationCount:Individuals:GrowthRate:Sterile  + 
#   MutationFrequency:Individuals:GrowthRate:Sterile )

# model.all.added <-glm(survivalRate~MutationCount*(I(MutationFrequency^2)+MutationFrequency)*(I(Individuals^2)+Individuals)*GrowthRate*Sterile,family=binomial(),data=groupedData)


library(ggplot2)



dataset <-"./data/out_LoadEquivalenceParameterSpace_LIFPNG_1753296009.csv"
title <- "linked"

groupedGraphData <- read.csv(dataset) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, extinctionRate, count) %>% 
  summarise()

# L_c
groupedGraphData$MeanIndividualMuts <- 2*groupedGraphData$MutationFrequency *groupedGraphData$MutationCount
# L_e
groupedGraphData$FilterRate <- 1-(1-groupedGraphData$MutationFrequency^2)^groupedGraphData$MutationCount

# gr<-3
# i<-25

# focusedDataByLoad <- groupedGraphData %>%
# filter(GrowthRate == gr) %>% 
# filter(Individuals == i)
colors <- c(
  '#a9db47',
  '#d9bd1e',
  '#d9ab13',
  '#d7990b',
  '#d38804',
  '#cd7700',
  '#c76500',
  '#c05300',
  '#b93f00',
  '#b02900',
  '#a80000')
  

dotpointscale = 0.75/1.35
# svglite("figures/figure_4_old.svg", width = 13*dotpointscale, height = 10*dotpointscale)
# par(mar=c(5,4,2,2)+0.1)
# ggplot(focusedDataByLoad,aes(x=MutationFrequency,y=MutationCount),) +
# xlab(expression("Deleterious recessive frequency (" * q * ")")) +
# ylab(expression("Deleterious loci count (" * l * ")")) +
# # ggtitle(sprintf("%s gr=%s i=%s",title,gr,i)) +
# geom_tile(aes(fill=extinctionRate)) +
# # scale_x_continuous(breaks = seq(0, max(focusedDataByLoad$Individuals), 25),expand = c(0,0)) +
# # scale_y_continuous(breaks = seq(0, 150, 5) ,expand = c(0,0)) +
# scale_fill_gradientn(name = "Extinction rate", colours = colors, breaks=c(0,0.5,1),limits=c(0,1)) +
# theme_classic() +
# labs(fill = "Generations") +
# theme(legend.title.align=0.5)
# dev.off()


svglite("figures/figure_4.svg", width = 13*dotpointscale, height = 18*dotpointscale)
par(mar=c(5,4,2,2)+0.1)

ggplot(groupedGraphData,aes(x=MutationFrequency,y=MutationCount),) +
  # ggtitle(sprintf("%s gr=%s i=%s",title,gr,i)) +
  geom_tile(aes(fill=extinctionRate)) +
  theme(panel.background = element_blank()) +
  scale_fill_gradientn(name = "Extinction rate", colours = colors, guide="none",limits=c(0,1)) +
  xlab(expression("Deleterious recessive frequency (" * q * ")")) +
  ylab(expression("Deleterious loci count (" * l * ")")) +
  theme_classic() +
  labs(fill = "Generations") +
  theme(legend.title.align=0.5) +
  facet_grid(vars(Individuals), vars(GrowthRate))
dev.off()

# linear model stats

for(gr in c(3,5)) {
  for(i in c(12,25,50,100)) {
    focusedDataByLoad <- groupedGraphData %>%
      filter(GrowthRate == gr) %>% 
      filter(Individuals == i)

    diff <- 0.475
    middlesubset <- filter(focusedDataByLoad,extinctionRate> 0.5-diff,extinctionRate< 0.5+diff)

    model<-glm(extinctionRate~FilterRate+MeanIndividualMuts,family = binomial(),data=middlesubset)
    print("-----------")
    print(gr)
    print(i)
    print(summary(model))

  }
}


# linear model AIC investivgation 
#  note that in some cases step() results in a different model, however in all cases the suggested model has a higher AIC than the above
for(gr in c(3,5)) {
  for(i in c(12,25,50,100)) {
    focusedDataByLoad <- groupedGraphData %>%
      filter(GrowthRate == gr) %>% 
      filter(Individuals == i)

    diff <- 0.475
    middlesubset <- filter(focusedDataByLoad,extinctionRate> 0.5-diff,extinctionRate< 0.5+diff)

    model<-glm(extinctionRate~FilterRate+MeanIndividualMuts+GrowthRate+Individuals+MutationFrequency+MutationCount,family = binomial(),data=middlesubset)
    step(model,direction='both')
    print(gr)
    print("-----------")
    print(i)
  }
}