library(dplyr)
# install.packages("svglite")
library(svglite)
# # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2 graphs
mc <- 200
mf <- 0.04
i <- 25
gr1 <- 3
gr2 <- 3.5
gr3 <- 4

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
dev.off()




# Analysis not included in the paper, compares loads that also affect both sexes vs females and sterile vs lethal of each.
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

