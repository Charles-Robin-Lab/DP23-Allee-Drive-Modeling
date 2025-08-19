library(dplyr)
# install.packages("svglite")
library(svglite)
# # groupedData <- read.csv("./data/out_51093850.csv") %>%
# #   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
# #   mutate(count = n()) %>%
# #   mutate(survivalRate = sum(Result != "EXTINCT") / count) %>%
# #   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
# #   summarise()
# # autosomalData <- groupedData[groupedData$Xlinked==0,]
# groupedGraphData <- read.csv("./data/out_GraphSlicesLoadTypeCompare_LIFP_1697426721.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, extinctionRate, count) %>% 
#   summarise()
# autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]
# # autosomalGraphData <- filter(autosomalGraphData,MutationFrequency!=0.5)

# # dataSlice2d.MutationFrequency <- filter(autosomalData, MutationCount==61, Individuals == 22, GrowthRate == 2, Sterile==1)
# # dataSlice2d.MutationCount <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==1)
# # dataSlice2d.Lethal <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==0)
# # dataSlice2d.Individuals <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, GrowthRate == 2, Sterile==1)
# # dataSlice2d.GrowthRate <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, Individuals == 22, Sterile==1)


# groupedGraphData <- read.csv("./data/out_GraphSlices100_LIFP_1704696372.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate, extinctionRate, count) %>% 
#   summarise()
# autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]

# gr <- 3
# mf <- 0.1
# dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == gr, Sterile==1)
# dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == 25, GrowthRate == gr, Sterile==1)
# dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==100, GrowthRate == gr, Sterile==1)

# # MutationCount==60, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4.0, Sterile==1

# dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==60, Individuals == 20, GrowthRate == 4, Sterile==0)
# dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4, Sterile==0)
# # dataSlice2d.Lethal <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 4, Sterile==0)
# dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, GrowthRate == 4, Sterile==0)
# dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, Individuals == 20, Sterile==0)



# # MutationCount==100, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3.0, Sterile==1

# dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == 3,Xlinked==0, FemaleOnlyEffect==1, Sterile==1)
# dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3, Sterile==1)
# dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, GrowthRate == 3, Sterile==1)
# dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, Individuals == 25, Sterile==1)

# dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==100, Individuals == 25, GrowthRate == 3, Xlinked==0, FemaleOnlyEffect==1, Sterile==0)
# dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.05, Individuals == 25, GrowthRate == 3, Sterile==0)
# dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, GrowthRate == 3, Sterile==0)
# dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.05, MutationCount==100, Individuals == 25, Sterile==0)
# # basic plots

# plot(extinctionRate~MutationFrequency,data=autosomalGraphData,ylim=c(0,1))

# plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency,ylim=c(0,1),ylab="Extinction rate",xlab=expression("Deleterious recessive frequency (" * q * ")"))
# plot(extinctionRate~MutationCount,data=dataSlice2d.MutationCount,ylim=c(0,1),ylab="Extinction rate",xlab=expression("Deleterious loci count (" * l * ")"))
# plot(extinctionRate~Individuals,data=dataSlice2d.Individuals,ylim=c(0,1),ylab="Extinction rate",xlab="Founding population size (N(0))")
# plot(extinctionRate~GrowthRate,data=dataSlice2d.GrowthRate,ylim=c(0,1),ylab="Extinction rate",xlab=expression("Female reproductive output (" * b[f] * ")"))

# # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2 graphs
mc <- 200
mf <- 0.04
i <- 25
gr1 <- 3
gr2 <- 3.5
gr3 <- 4

# groupedGraphData <- read.csv("./data/out_GraphSlices80.3_LIFP_1735381152.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
#   summarise()
# autosomalGraphData <- filter(groupedGraphData,Xlinked==0,Sterile==0)
# dataSlice2d.MutationFrequency1 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr1)
# dataSlice2d.MutationCount1 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr1)
# dataSlice2d.Individuals1 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr1)
# dataSlice2d.GrowthRate1 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)

# groupedGraphData <- read.csv("./data/out_GraphSlices80.3.5_LIFP_1707968395.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
#   summarise()
# autosomalGraphData <- filter(groupedGraphData,Xlinked==0,Sterile==0)
# dataSlice2d.MutationFrequency2 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr2)
# dataSlice2d.MutationCount2 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr2)
# dataSlice2d.Individuals2 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr2)
# dataSlice2d.GrowthRate2 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)

# groupedGraphData <- read.csv("./data/out_GraphSlices80_LIFP_1707968384.csv") %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
#   summarise()
# autosomalGraphData <- filter(groupedGraphData,Xlinked==0,Sterile==0)
# dataSlice2d.MutationCount3 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr3)
# dataSlice2d.MutationFrequency3 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr3)
# dataSlice2d.Individuals3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr3)
# dataSlice2d.GrowthRate3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)

# groupedGraphData <- read.csv("./data/out_GraphSlicesAllStandardLinkage_LIFPIC_1747151288.csv") %>%
# groupedGraphData <- read.csv("./data/out_GraphSlicesAllStandardLinkage400_LIFPIC_1747236301.csv") %>%
groupedGraphData <- read.csv("./data/out_GraphSlicesAll_LIFPNG_1753085746.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, PostCompetitionMutationTiming) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, PostCompetitionMutationTiming, extinctionRate, count) %>% 
  summarise()
autosomalGraphData <- filter(groupedGraphData,Xlinked==0,Sterile==0 |Sterile=="F",PostCompetitionMutationTiming == 0) 
dataSlice2d.MutationFrequency1 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr1,MutationFrequency<=0.07)
dataSlice2d.MutationCount1 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr1)
dataSlice2d.Individuals1 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr1)
dataSlice2d.GrowthRate1 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)
dataSlice2d.MutationFrequency2 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr2,MutationFrequency<=0.07)
dataSlice2d.MutationCount2 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr2)
dataSlice2d.Individuals2 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr2)
dataSlice2d.GrowthRate2 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)
dataSlice2d.MutationCount3 <- filter(autosomalGraphData, MutationFrequency==mf, Individuals == i, GrowthRate == gr3)
dataSlice2d.MutationFrequency3 <- filter(autosomalGraphData, MutationCount==mc, Individuals == i, GrowthRate == gr3,MutationFrequency<=0.07)
dataSlice2d.Individuals3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, GrowthRate == gr3)
dataSlice2d.GrowthRate3 <- filter(autosomalGraphData, MutationFrequency==mf, MutationCount==mc, Individuals == i)



svglite("figures/figure_2A.svg", width = 5.8, height = 4.35)
plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate1,GrowthRate!=gr1,GrowthRate!=gr2,GrowthRate!=gr3,2<GrowthRate*(1-mf^2)^mc),ylim=c(0,1),xlim=c(2,6),ylab="Extinction rate",xlab=expression("Female reproductive output (" * b[f] * ")"),pch=19)
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate1,GrowthRate!=gr1,GrowthRate!=gr2,GrowthRate!=gr3,2>=GrowthRate*(1-mf^2)^mc),ylim=c(0,1),xlim=c(2,6),axes=FALSE,pch=1,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate1,GrowthRate==gr1),ylim=c(0,1),xlim=c(2,6),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate2,GrowthRate==gr2),ylim=c(0,1),xlim=c(2,6),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate3,GrowthRate==gr3),ylim=c(0,1),xlim=c(2,6),axes=FALSE,col="#68c24f",pch=19,ylab="",xlab="")
legend("bottomleft", inset = c(0.05, 0.05),  
  title=expression(b[f]),
  legend = c(3,3.5,4),
  pch=16,
  col=c("#81651d","#859225","#68c24f"))
dev.off()

svglite("figures/figure_2B.svg", width = 5.8, height = 4.35)
# plot(extinctionRate~MutationFrequency,data=dataSlice2d.MutationFrequency4,ylim=c(0,1),xlim=c(0,0.15),ylab="Extinction rate",xlab=expression("Deleterious recessive frequency (" * q * ")"),col="#ff0000",pch=19)
# par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency3,2>=gr3*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),ylab="Extinction rate",xlab=expression("Deleterious recessive frequency (" * q * ")"),col="#68c24f",pch=1)
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency2,2>=gr2*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),axes=FALSE,col="#859225",pch=1,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency1,2>=gr1*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),axes=FALSE,col="#81651d",pch=1,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency3,2<gr3*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),axes=FALSE,col="#68c24f",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency2,2<gr2*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationFrequency,data=filter(dataSlice2d.MutationFrequency1,2<gr1*(1-MutationFrequency^2)^mc),ylim=c(0,1),xlim=c(0,0.07),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
abline(v=mf,lty=2)
legend("bottomright", inset = c(0.05, 0.05),  
  title=expression(b[f]),
  legend = c(3,3.5,4),
  pch=16,
  col=c("#81651d","#859225","#68c24f"))
dev.off()

svglite("figures/figure_2C.svg", width = 5.8, height = 4.35)
# plot(extinctionRate~Individuals,data=dataSlice2d.Individuals4,ylim=c(0,1),xlim=c(0,0.15),ylab="Extinction rate",xlab=expression("Deleterious recessive frequency (" * q * ")"),col="#68c24f",pch=19)
# par(new=TRUE)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals3,ylim=c(0,1),xlim=c(0,120),ylab="Extinction rate",xlab="Founding population size (N(0))",col="#68c24f",pch=19)
par(new=TRUE)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals2,ylim=c(0,1),xlim=c(0,120),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~Individuals,data=dataSlice2d.Individuals1,ylim=c(0,1),xlim=c(0,120),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
abline(v=i,lty=2)
dev.off()

svglite("figures/figure_2D.svg", width = 5.8, height = 4.35)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount1,2>=gr1*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),axes=FALSE,col="#81651d",pch=1,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount2,2>=gr2*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),axes=FALSE,col="#859225",pch=1,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount3,2>=gr3*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),ylab="Extinction rate",xlab=expression("Deleterious loci count (" * l * ")"),col="#68c24f",pch=1)
par(new=TRUE)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount1,2<gr1*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),axes=FALSE,col="#81651d",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount2,2<gr2*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),axes=FALSE,col="#859225",pch=19,ylab="",xlab="")
par(new=TRUE)
plot(extinctionRate~MutationCount,data=filter(dataSlice2d.MutationCount3,2<gr3*(1-mf^2)^MutationCount),ylim=c(0,1),xlim=c(0,400),axes=FALSE,ylab="",xlab="",col="#68c24f",pch=19)
abline(v=mc,lty=2)
legend("bottomright", inset = c(0.05, 0.05),  
  title=expression(b[f]),
  legend = c(3,3.5,4),
  pch=16,
  col=c("#81651d","#859225","#68c24f"))
dev.off()




# Analysis not included in the paper, compares loads that also affect males.
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
