# entire sim fitness trends
library(dplyr)
library(svglite)
FitnessData <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
FitnessData90_4 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_66x4.csv")
FitnessData85_4 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
FitnessData90_2 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_264x2.csv")
FitnessData85_2 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_407x2.csv")
carryingCapacity <- 1000;
FitnessData$percentOfFitness <- FitnessData$Fitness / ( carryingCapacity * FitnessData$GrowthRate / (carryingCapacity + FitnessData$GrowthRate * FitnessData$Females - FitnessData$Individuals))
FitnessData90_4$percentOfFitness <- FitnessData90_4$Fitness / ( carryingCapacity * FitnessData90_4$GrowthRate / (carryingCapacity + FitnessData90_4$GrowthRate * FitnessData90_4$Females - FitnessData90_4$Individuals))
FitnessData85_4$percentOfFitness <- FitnessData85_4$Fitness / ( carryingCapacity * FitnessData85_4$GrowthRate / (carryingCapacity + FitnessData85_4$GrowthRate * FitnessData85_4$Females - FitnessData85_4$Individuals))
FitnessData90_2$percentOfFitness <- FitnessData90_2$Fitness / ( carryingCapacity * FitnessData90_2$GrowthRate / (carryingCapacity + FitnessData90_2$GrowthRate * FitnessData90_2$Females - FitnessData90_2$Individuals))
FitnessData85_2$percentOfFitness <- FitnessData85_2$Fitness / ( carryingCapacity * FitnessData85_2$GrowthRate / (carryingCapacity + FitnessData85_2$GrowthRate * FitnessData85_2$Females - FitnessData85_2$Individuals))

any(filter(FitnessData90_4,Time==2)$count!=1)
# FitnessData <- read.csv("./data/out_GraphSlices80.3_LIFP_1735381152.csv") %>%
#   group_by(GrowthRate, Time) %>% 
#   mutate(count = n()) %>%
#   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
#   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
#   summarise()

firstGeneration <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv")
params <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv_params.csv")
focusedSeeds <- filter(params,Sterile==0,FemaleOnlyEffect==1)$Seed
FitnessData <- filter(firstGeneration, Seed %in% focusedSeeds)
# FitnessData <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
FitnessData$percentOfFitness <- FitnessData$Fitness / ( carryingCapacity * FitnessData$GrowthRate / (carryingCapacity + FitnessData$GrowthRate * FitnessData$Females - FitnessData$Individuals))

lastSize <- FitnessData %>%
  group_by(Seed) %>% 
  mutate(lastGeneration = max(Time))  %>%
  filter(lastGeneration==Time) %>%
  group_by(Seed, Time, GrowthRate, lastGeneration, Individuals, Females) %>% 
  summarise()
hist(lastSize$Individuals,nclass=100)  

sum(lastSize$Individuals < 500)/(sum(lastSize$Individuals > 500) + sum(lastSize$Individuals <= 500))


survivedSims <- FitnessData %>%
  filter(Seed %in% lastSize$Seed[lastSize$Individuals > 500])
deadSims <- FitnessData %>%
  filter(Seed %in% lastSize$Seed[lastSize$Individuals <= 500])

hist(lastSize$lastGeneration[lastSize$Individuals > 500],breaks = seq(0, 150, 2))
hist(lastSize$lastGeneration[lastSize$Individuals <= 500],breaks = seq(0, 60, 1))

lowGrowthIndividuals <- survivedSims %>%
  filter(GrowthRate == 3) %>% 
  filter(Time == 3) %>% 
  group_by(Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
midGrowthIndividuals <- survivedSims %>%
  filter(GrowthRate == 3.5) %>% 
  group_by(Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
highGrowthIndividuals <- survivedSims %>%
  filter(GrowthRate == 4) %>% 
  group_by(Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()

lowGrowthTime <- deadSims %>%
  filter(GrowthRate == 3) %>% 
  group_by(Time) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
midGrowthTime <- survivedSims %>%
  filter(GrowthRate == 3.5) %>% 
  group_by(Time) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
highGrowthTime <- survivedSims %>%
  filter(GrowthRate == 4) %>% 
  group_by(Time) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
  summarise()

lowGrowthTimeAndInds <- survivedSims %>%
  filter(GrowthRate == 3) %>% 
  filter(percentOfFitness != 0) %>% 
  group_by(Time, Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
midGrowthTimeAndInds <- survivedSims %>%
  filter(GrowthRate == 3.5) %>% 
  group_by(Time, Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()
highGrowthTimeAndInds <- deadSims %>%
  filter(GrowthRate == 4) %>% 
  group_by(Time, Individuals) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
  summarise()

plot(meanPercentOfFitness~Time,data=lowGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
plot(meanPercentOfFitness~Time,data=midGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
plot(meanPercentOfFitness~Time,data=highGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
plot(meanPercentOfFitness~Individuals,data=lowGrowthIndividuals)
plot(meanPercentOfFitness~Individuals,data=midGrowthIndividuals)
plot(meanPercentOfFitness~Individuals,data=highGrowthIndividuals)

library(lattice)
wireframe(meanPercentOfFitness ~ Individuals * Time, data=lowGrowthTimeAndInds)
wireframe(meanPercentOfFitness ~ Individuals * Time, data=midGrowthTimeAndInds)
wireframe(meanPercentOfFitness ~ Individuals * Time, data=highGrowthTimeAndInds)

library(ggplot2)
library(scales)

colors <- c("#ff0000", 
            "#ff0000",
            "#ff0000", 
            "#ff0000",  
            "#ff0000", 
            "#00ffff", 
            "#ffff00", 
            "#0000ff", 
            "#00ff00", 
            "#000000")

ggplot(lowGrowthTimeAndInds,aes(x=Individuals,y=Time),) +
  ggtitle("survived") +
  geom_tile(aes(fill=meanPercentOfFitness)) +
  scale_x_continuous(breaks = seq(0, max(lowGrowthTimeAndInds$Individuals), 25),expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 150, 5) ,expand = c(0,0)) +
  scale_fill_gradientn(name = "", colours = colors ,guide="none") +
  theme_classic() +
  labs(fill = "Generations") +
  theme(legend.title.align=0.5)



# 
library(dplyr)
firstGeneration <- read.csv("./data/out_FitnessInFirstGen_LIFP_1744373846.csv")

carryingCapacity <- 1000;
firstGeneration$percentOfFitness <- firstGeneration$Fitness / ( carryingCapacity * firstGeneration$GrowthRate / (carryingCapacity + firstGeneration$GrowthRate * firstGeneration$Females - firstGeneration$Individuals))

lowGrowthfirstGenerationIF <- firstGeneration %>%
  filter(GrowthRate == 3) %>% 
  group_by(Individuals, Females) %>% 
  mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
  group_by(Individuals, Females, GrowthRate, meanPercentOfFitness) %>% 
  summarise()

# lowGrowthfirstGenerationI <- firstGeneration %>%
#   filter(GrowthRate == 3) %>% 
#   group_by(Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()

# lowGrowthfirstGenerationF <- firstGeneration %>%
#   filter(GrowthRate == 3) %>% 
#   group_by(Females) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Females, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()


library(ggplot2)
library(scales)

breaks <- c(0,0.5,0.8,0.9,0.95,1.0)
colors <- c("#ff0000", "#ff00ff", "#00ff00", "#00ffff", "#0000ff", "#000000")

# Create a named vector of colors for the scale
color_scale <- setNames(colors, formatC(breaks, format = "g"))


ggplot(lowGrowthfirstGenerationIF,aes(x=Individuals,y=Females),) +
  geom_tile(aes(fill=meanPercentOfFitness)) +
  scale_x_continuous(breaks = seq(0, max(firstGeneration$Individuals), 5),expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, max(firstGeneration$Females), 5) ,expand = c(0,0)) +
  scale_fill_gradientn(name = "", colours = color_scale, breaks = breaks ,guide="none") +
  theme_classic() +
  labs(fill = "Generations") +
  theme(legend.title.align=0.5)


plot(meanPercentOfFitness~Females,data=lowGrowthfirstGenerationF,ylim=c(0.5,1),xlim=c(0,50))


plot(meanPercentOfFitness~Individuals,data=lowGrowthfirstGenerationI,ylim=c(0.5,1),xlim=c(0,65))



library(dplyr)
library(svglite)
firstGeneration <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1750493264.csv")
params <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1750493264.csv_params.csv")
focusedSeeds <- filter(params,Sterile==1,FemaleOnlyEffect==1)$Seed
focusedFirstGeneration <- filter(firstGeneration, Seed %in% focusedSeeds)
# firstGeneration <- FitnessData90_2
preSorted <- focusedFirstGeneration %>% 
    group_by(Seed) %>% 
    arrange(Time)

svglite("figures/figure_2.53s.svg", width = 8, height = 6)
colors <- c("#68c24f","#859225", "#81651d","#811d1d","#591414","#a000f0", "#000000")
points <- 1+c(1,3,5,10,15,25,40)
# points <- 1+c(1,2,3,4,5,7,10)
# points <- 1+c(1,10,11,12,13,14,15)
# points <- 1+c(1,16,20,24,28,32,36)

for (j in 1:length(points)) {
  i <- points[j]
  nthGeneration <- filter(preSorted,Time<=i)

  survivedSims <- nthGeneration %>% 
    filter(n()==i-1) %>%
    summarise(
      finalPercentOfFitness = last(Fitness),
      startSize = first(Individuals),
      startFemales = first(Females),
      secondLastFemales = nth(Females,i-2),
      secondLastSize = nth(Individuals,i-2),
      finalFemales = nth(Females,i-1),
      finalSize = nth(Individuals,i-1),
      GrowthRate
    ) %>%
    filter(finalSize!=finalFemales) %>%
    filter(0!=finalFemales)

  lowGrowthSecondGenerationI <- survivedSims %>%
    filter(GrowthRate == 3) %>% 
    group_by(startSize) %>% 
    mutate(meanFinalPercentOfFitness = mean(finalPercentOfFitness))  %>%
    group_by(startSize, meanFinalPercentOfFitness) %>% 
    summarise()
  # lowGrowthSecondGenerationI <- survivedSims %>%
  #   filter(GrowthRate == 3) %>% 
  #   group_by(startSize) %>% 
  #   mutate(meanFinalPercentOfFitness = mean(finalPercentOfFitness))  %>%
  #   group_by(startSize, meanFinalPercentOfFitness) %>% 
  #   summarise()

  plot(meanFinalPercentOfFitness~startSize,data=lowGrowthSecondGenerationI,
    xlim=c(0,80),
    ylim=c(0.0,1.0),
    col=colors[j],
    pch=16,
    axes=j==1,
    type = "o",
    ylab=if (j==1) "Normalised mean fitness" else "",
    xlab=if (j==1) "Founding population size" else "")
  par(new=TRUE)

}
par(new=FALSE)
legend("bottomright", inset = c(0.05, 0.05),  
  title="Generations",
  legend = points-1,
  pch=16,
  col=colors)
dev.off()

params <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv_params.csv")
groupedGraphData <- params %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, RecombinationRate, extinctionRate, count) %>% 
  summarise()
focusedSims <- filter(groupedGraphData,Sterile==1,FemaleOnlyEffect==1)
focusedSims2 <- filter(groupedGraphData,Sterile==0,FemaleOnlyEffect==1)
plot(extinctionRate~Individuals,data=focusedSims,ylim=c(0,1),ylab="Extinction probability",xlab="Founding population size")
par(new=TRUE)
plot(extinctionRate~Individuals,data=focusedSims2,ylim=c(0,1),ylab="Extinction probability",xlab="Founding population size",col="#ff0000")

plot(extinctionRate~GrowthRate,data=filter(dataSlice2d.GrowthRate3,GrowthRate==4),ylim=c(0,1),xlim=c(2,6),axes=FALSE,pch=19,ylab="",xlab="")



