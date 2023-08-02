library(dplyr)

groupedData <- read.csv("../data/out.csv") %>%
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

