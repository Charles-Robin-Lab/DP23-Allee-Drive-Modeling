library(dplyr)
library(svglite)

# groupedData <- read.csv("./data/out_51093850.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
#   mutate(count = n()) %>%
#   mutate(survivalRate = sum(Result != "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
#   summarise()
# autosomalData <- groupedData[groupedData$Xlinked==0,]
groupedGraphData <- read.csv("./data/out_GraphSlicesLoadTypeCompare_LIFP_1697426721.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]
# autosomalGraphData <- filter(autosomalGraphData,MutationFrequency!=0.5)

# dataSlice2d.MutationFrequency <- filter(autosomalData, MutationCount==61, Individuals == 22, GrowthRate == 2, Sterile==1)
# dataSlice2d.MutationCount <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==1)
# dataSlice2d.Lethal <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==0)
# dataSlice2d.Individuals <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, GrowthRate == 2, Sterile==1)
# dataSlice2d.GrowthRate <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, Individuals == 22, Sterile==1)


groupedGraphData <- read.csv("./data/out_GraphSlices100_LIFP_1704696372.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]

gr <- 3
mf <- 0.1
dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == gr, Sterile==1)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == 25, GrowthRate == gr, Sterile==1)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==100, GrowthRate == gr, Sterile==1)

# MutationCount==60, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4.0, Sterile==1

dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==60, Individuals == 20, GrowthRate == 4, Sterile==0)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4, Sterile==0)
# dataSlice2d.Lethal <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4, Sterile==0)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, GrowthRate == 4, Sterile==0)
dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, Individuals == 20, Sterile==0)



# MutationCount==100, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3.0, Sterile==1

dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == 3,Xlinked==0, FemaleOnlyEffect==1, Sterile==1)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3, Sterile==1)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, GrowthRate == 3, Sterile==1)
dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, Individuals == 25, Sterile==1)

dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == 3, Xlinked==0, FemaleOnlyEffect==1, Sterile==0)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3, Sterile==0)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, GrowthRate == 3, Sterile==0)
dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, Individuals == 25, Sterile==0)
# basic plots

plot(extinctionRate~MutationFrequency,data=autosomalGraphData,ylim=c(0,1))

plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious recessive frequency")
plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious loci count")
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals,ylim=c(0,1),ylab="Extinction probability",xlab="Founding population size")
plot(extinctionRate~GrowthRate,data=dataSlice2d.GrowthRate,ylim=c(0,1),ylab="Extinction probability",xlab="Female reproductive output ")

# # # # # # # # # # # # # # # # # # # # # # # # 
# find the cool graphs
mc <- 80
mf <- 0.075
i <- 25
gr1 <- 3
gr2 <- 3.5
gr3 <- 4

groupedGraphData <- read.csv("./data/out_GraphSlices80,3_LIFP_1707892121.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]
dataSlice2d.MutationFrequency1 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr1, Sterile==0)
dataSlice2d.MutationCount1 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr1, Sterile==0)
dataSlice2d.Individuals1 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr1, Sterile==0)

groupedGraphData <- read.csv("./data/out_GraphSlices80.3.5_LIFP_1707968395.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]
dataSlice2d.MutationFrequency2 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr2, Sterile==0)
dataSlice2d.MutationCount2 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr2, Sterile==0)
dataSlice2d.Individuals2 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr2, Sterile==0)

groupedGraphData <- read.csv("./data/out_GraphSlices80_LIFP_1707968384.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]
dataSlice2d.MutationCount3 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr3, Sterile==0)
dataSlice2d.MutationFrequency3 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr3, Sterile==0)
dataSlice2d.Individuals3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr3, Sterile==0)
dataSlice2d.GrowthRate3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i, Sterile==0)

# dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==81, Individuals == 12, GrowthRate == 3.2, Xlinked==0, Sterile==0)
# dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.0, Individuals == 12, GrowthRate == 3.2, Xlinked==0, Sterile==0)
# dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.06, MutationCount==81, GrowthRate == 3.2, Xlinked==0, Sterile==0)
# dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.06, MutationCount==81, Individuals == 12, Xlinked==0, Sterile==0)



svglite("figures/figure_2A.svg", width = 5.8, height = 4.35)
plot(extinctionRate~GrowthRate,data=dataSlice2d.GrowthRate3,ylim=c(0,1),ylab="Extinction probability",xlab="Female reproductive output ",pch=19)
dev.off()

svglite("figures/figure_2B.svg", width = 5.8, height = 4.35)
plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency3,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious recessive frequency",col="#68c24f",pch=19)
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency2,ylim=c(0,1),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency1,ylim=c(0,1),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
dev.off()

svglite("figures/figure_2C.svg", width = 5.8, height = 4.35)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals3,ylim=c(0,1),ylab="Extinction probability",xlab="Founding population size",col="#68c24f",pch=19)
par(new=TRUE)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals2,ylim=c(0,1),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals1,ylim=c(0,1),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
dev.off()

svglite("figures/figure_2D.svg", width = 5.8, height = 4.35)
plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount3,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious loci count",col="#68c24f",pch=19)
par(new=TRUE)
plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount2,ylim=c(0,1),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount1,ylim=c(0,1),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
dev.off()

plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency2,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious recessive frequency",pch=19)
plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount2,ylim=c(0,1),ylab="Extinction probability",xlab="Deleterious loci count",pch=19)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals2,ylim=c(0,1),ylab="Extinction probability",xlab="Founding population size",pch=19)







##### Load comparison ####
groupedGraphData <- read.csv("./data/out_GraphSlicesLoadTypeCompare_LIFP_1697426721.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result != "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, survivalRate, count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]

dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==60, Individuals == 20, GrowthRate == 4)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4)
# dataSlice2d.Lethal <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, GrowthRate == 4)
dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, Individuals == 20)


plot(survivalRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency,Sterile==0,FemaleOnlyEffect==1),col="black",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency,Sterile==0,FemaleOnlyEffect==0),col="blue",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency,Sterile==1,FemaleOnlyEffect==1),col="red",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency,Sterile==1,FemaleOnlyEffect==0),col="purple",ylim=c(0,1))

plot(survivalRate~MutationCount,data=filter(dataSlice2d.MutationCount,Sterile==0,FemaleOnlyEffect==1),col="black",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationCount,data=filter(dataSlice2d.MutationCount,Sterile==0,FemaleOnlyEffect==0),col="blue",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationCount,data=filter(dataSlice2d.MutationCount,Sterile==1,FemaleOnlyEffect==1),col="red",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~MutationCount,data=filter(dataSlice2d.MutationCount,Sterile==1,FemaleOnlyEffect==0),col="purple",ylim=c(0,1))

plot(survivalRate~Individuals,data=filter(dataSlice2d.Individuals,Sterile==0,FemaleOnlyEffect==1),col="black",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~Individuals,data=filter(dataSlice2d.Individuals,Sterile==0,FemaleOnlyEffect==0),col="blue",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~Individuals,data=filter(dataSlice2d.Individuals,Sterile==1,FemaleOnlyEffect==1),col="red",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~Individuals,data=filter(dataSlice2d.Individuals,Sterile==1,FemaleOnlyEffect==0),col="purple",ylim=c(0,1))

plot(survivalRate~GrowthRate,data=filter(dataSlice2d.GrowthRate,Sterile==0,FemaleOnlyEffect==1),col="black",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~GrowthRate,data=filter(dataSlice2d.GrowthRate,Sterile==0,FemaleOnlyEffect==0),col="blue",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~GrowthRate,data=filter(dataSlice2d.GrowthRate,Sterile==1,FemaleOnlyEffect==1),col="red",ylim=c(0,1))
par(new=TRUE)
plot(survivalRate~GrowthRate,data=filter(dataSlice2d.GrowthRate,Sterile==1,FemaleOnlyEffect==0),col="purple",ylim=c(0,1))


# 3d ideas
library(ggplot2)
library(scales)


dataSlice3d.countfreq <- filter(groupedData, Individuals == 32, GrowthRate == 2.2, Sterile==1, Xlinked==1)

dataSlice3d.countInd <- filter(groupedData, MutationFrequency == 0.11, GrowthRate == 2.2, Sterile==1, Xlinked==1)

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
