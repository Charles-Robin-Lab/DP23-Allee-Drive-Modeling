
library(dplyr)

groupedData <- read.csv("./out.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()



XlessgroupedData <- read.csv("./newReproductionGraphs.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, survivalRate,count) %>% 
  summarise()


recombinationComparisonData <- read.csv("./recombinationComparison.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(anySurvivalRate = sum(Result == "SURVIVED" | Result == "LOADED_SURVIVAL") / count) %>%
  mutate(loadSurvivalRate = sum(Result == "LOADED_SURVIVAL") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, RecombinationRate, anySurvivalRate, loadSurvivalRate,count) %>% 
  summarise()

plot(anySurvivalRate~RecombinationRate,data=recombinationComparisonData,ylim=c(0,1),log="x")
par(new=TRUE)
plot(loadSurvivalRate~RecombinationRate,data=recombinationComparisonData,pch = 19,ylim=c(0,1),log="x")
par(new=TRUE)
plot(anySurvivalRate-loadSurvivalRate~RecombinationRate,data=recombinationComparisonData,pch = 2,ylim=c(0,1),log="x")
  

xtestData <- read.csv("./out.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  filter(survivalRate<0.975 && survivalRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()

xlinkedDiffs <- xtestData[xtestData$Xlinked==1,]$survivalRate-xtestData[xtestData$Xlinked==0,]$survivalRate
hist(xlinkedDiffs,breaks = c(-0.3,-0.2,-0.15,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.15,0.2,0.3))

steriletestData <- read.csv("./out.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  filter(survivalRate<0.975 && survivalRate>0.025) %>%
  ungroup()  %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()
sterileDiffs <- steriletestData[steriletestData$Sterile==1,]$survivalRate-steriletestData[steriletestData$Sterile==0,]$survivalRate
hist(sterileDiffs,breaks = c(-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.0375,0.0375,0.1,0.2,0.3))
boxplot(sterileDiffs,xlinkedDiffs)


XlessgroupedData <- XlessgroupedData %>%
  group_by(MutationFrequency, MutationCount, Individuals) %>% 
  mutate(count = n()) %>%
  mutate(groupedSurvivalRate = sum(survivalRate)/count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, groupedSurvivalRate) %>% 
  summarise()



# Extra stats
groupedData$logisticSurvivalRate <- log(groupedData$survivalRate/(1-groupedData$survivalRate))
groupedData$sQIndividuals <- groupedData$Individuals^2
groupedData$sQMutationFrequency <- groupedData$MutationFrequency^2
groupedData$eIndividuals <- exp(groupedData$Individuals)
groupedData$individualChance <- 0.5*(1-(1-groupedData$MutationFrequency^2)^groupedData$MutationCount)

# Compare X-linked and not
# dataSlice2d.mutationCount <- filter(groupedData, MutationFrequency==0.16, Individuals == 12, GrowthRate == 22, Sterile==1)
# dataSlice2d.mutationCount <- filter(groupedData, MutationFrequency==0.21, Individuals == 102, GrowthRate == 14, Xlinked==1)
# dataSlice2d.mutationCount <- filter(XlessgroupedData, MutationFrequency==0.11, Individuals == 42, GrowthRate == 2, Sterile==1)
dataSlice2d.MutationFrequency <- filter(XlessgroupedData, MutationCount==61, Individuals == 22, GrowthRate == 2, Sterile==1)
dataSlice2d.MutationCount <- filter(XlessgroupedData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==1)
dataSlice2d.Lethal <- filter(XlessgroupedData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==0)
dataSlice2d.Individuals <- filter(XlessgroupedData, MutationFrequency==0.11, MutationCount==61, GrowthRate == 2, Sterile==1)
dataSlice2d.GrowthRate <- filter(XlessgroupedData, MutationFrequency==0.11, MutationCount==61, Individuals == 22, Sterile==1)




# MutationCount==60, MutationFrequency==0.10, Individuals == 20, GrowthRate == 1, Sterile==1

dataSlice2d.MutationFrequency <- filter(XlessgroupedData, MutationCount==60, Individuals == 20, GrowthRate == 1, Sterile==1)
dataSlice2d.MutationCount <- filter(XlessgroupedData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 1, Sterile==1)
dataSlice2d.Lethal <- filter(XlessgroupedData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 2, Sterile==0)
dataSlice2d.Individuals <- filter(XlessgroupedData, MutationFrequency==0.10, MutationCount==60, GrowthRate == 1, Sterile==1)
dataSlice2d.GrowthRate <- filter(XlessgroupedData, MutationFrequency==0.10, MutationCount==60, Individuals == 20, Sterile==1)



plot(survivalRate~MutationFrequency,data=dataSlice2d.MutationFrequency,ylim=c(0,1))
plot(survivalRate~MutationCount,data=dataSlice2d.MutationCount,ylim=c(0,1))
plot(survivalRate~Individuals,data=dataSlice2d.Individuals,ylim=c(0,1))
plot(survivalRate~GrowthRate,data=dataSlice2d.GrowthRate,ylim=c(0,1))




library(ggplot2)
dataSlice3d.countfreq <- filter(groupedData, Individuals == 32, GrowthRate == 2, Sterile==1, Xlinked==1)

dataSlice3d.countInd <- filter(groupedData, MutationFrequency == 0.11, GrowthRate == 2, Sterile==1, Xlinked==1)

breaks <- c(22,24,28,36,52,84,120)
colors <- c("#fcfbc7", "#efb87c", "#e56a57" , "#be3371", "#721f7d", "#300c61", "#000000")

# Create a named vector of colors for the scale
color_scale <- setNames(colors, formatC(breaks, format = "g"))


ggplot(dataSlice3d.countfreq,aes(x=MutationFrequency,y=MutationCount),)+
  geom_tile(aes(fill=survivalRate))+
  scale_x_continuous(breaks = seq(0.01, 0.31, 0.05),expand = c(0,0) )+
  scale_y_continuous(breaks = seq(2, 102, 10) ,expand = c(0,0))+
  scale_fill_gradientn(name = "", colours = color_scale, breaks = breaks ,guide="none") +
  theme_classic()+
  labs(fill = "Generations")+
  theme(legend.title.align=0.5)
  
ggplot(dataSlice3d.countInd,aes(x=Individuals,y=MutationCount),)+
  geom_tile(aes(fill=survivalRate))+
  scale_x_continuous(breaks = seq(2, 102, 10),expand = c(0,0) )+
  scale_y_continuous(breaks = seq(2, 102, 10) ,expand = c(0,0))+
  scale_fill_gradientn(name = "", colours = color_scale, breaks = breaks ,guide="none") +
  theme_classic()+
  labs(fill = "Generations")+
  theme(legend.title.align=0.5)


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


dataLine <- filter(groupedData,survivalRate <= 0.90,survivalRate >= 0.10)

# time analysis
x = (0:(1500/25))*25
y= data.frame(w=0)
for (a in x) {
   y[nrow(y)+1,] = nrow(data[data$Time>a,])/nrow(data)
}
plot(x,y[2:62,])
# print(filter(data, Time>5000, Result == "EXTINCT"))