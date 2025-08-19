# # entire sim fitness trends
# library(dplyr)
# library(svglite)
# FitnessData <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
# FitnessData90_4 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_66x4.csv")
# FitnessData85_4 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
# FitnessData90_2 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_264x2.csv")
# FitnessData85_2 <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_407x2.csv")
# carryingCapacity <- 1000;
# FitnessData$percentOfFitness <- FitnessData$Fitness / ( carryingCapacity * FitnessData$GrowthRate / (carryingCapacity + FitnessData$GrowthRate * FitnessData$Females - FitnessData$Individuals))
# FitnessData90_4$percentOfFitness <- FitnessData90_4$Fitness / ( carryingCapacity * FitnessData90_4$GrowthRate / (carryingCapacity + FitnessData90_4$GrowthRate * FitnessData90_4$Females - FitnessData90_4$Individuals))
# FitnessData85_4$percentOfFitness <- FitnessData85_4$Fitness / ( carryingCapacity * FitnessData85_4$GrowthRate / (carryingCapacity + FitnessData85_4$GrowthRate * FitnessData85_4$Females - FitnessData85_4$Individuals))
# FitnessData90_2$percentOfFitness <- FitnessData90_2$Fitness / ( carryingCapacity * FitnessData90_2$GrowthRate / (carryingCapacity + FitnessData90_2$GrowthRate * FitnessData90_2$Females - FitnessData90_2$Individuals))
# FitnessData85_2$percentOfFitness <- FitnessData85_2$Fitness / ( carryingCapacity * FitnessData85_2$GrowthRate / (carryingCapacity + FitnessData85_2$GrowthRate * FitnessData85_2$Females - FitnessData85_2$Individuals))

# any(filter(FitnessData90_4,Time==2)$count!=1)
# # FitnessData <- read.csv("./data/out_GraphSlices80.3_LIFP_1735381152.csv") %>%
# #   group_by(GrowthRate, Time) %>% 
# #   mutate(count = n()) %>%
# #   mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
# #   group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, RecombinationRate, extinctionRate, count) %>% 
# #   summarise()

# firstGeneration <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv")
# params <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv_params.csv")
# focusedSeeds <- filter(params,Sterile==0,FemaleOnlyEffect==1)$Seed
# FitnessData <- filter(firstGeneration, Seed %in% focusedSeeds)
# # FitnessData <- read.csv("./data/out_FitnessGivenDifferentLoads30_LIFP_1745690157.csv_102x4.csv")
# FitnessData$percentOfFitness <- FitnessData$Fitness / ( carryingCapacity * FitnessData$GrowthRate / (carryingCapacity + FitnessData$GrowthRate * FitnessData$Females - FitnessData$Individuals))

# lastSize <- FitnessData %>%
#   group_by(Seed) %>% 
#   mutate(lastGeneration = max(Time))  %>%
#   filter(lastGeneration==Time) %>%
#   group_by(Seed, Time, GrowthRate, lastGeneration, Individuals, Females) %>% 
#   summarise()
# hist(lastSize$Individuals,nclass=100)  

# sum(lastSize$Individuals < 500)/(sum(lastSize$Individuals > 500) + sum(lastSize$Individuals <= 500))


# survivedSims <- FitnessData %>%
#   filter(Seed %in% lastSize$Seed[lastSize$Individuals > 500])
# deadSims <- FitnessData %>%
#   filter(Seed %in% lastSize$Seed[lastSize$Individuals <= 500])

# hist(lastSize$lastGeneration[lastSize$Individuals > 500],breaks = seq(0, 150, 2))
# hist(lastSize$lastGeneration[lastSize$Individuals <= 500],breaks = seq(0, 60, 1))

# lowGrowthIndividuals <- survivedSims %>%
#   filter(GrowthRate == 3) %>% 
#   filter(Time == 3) %>% 
#   group_by(Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# midGrowthIndividuals <- survivedSims %>%
#   filter(GrowthRate == 3.5) %>% 
#   group_by(Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# highGrowthIndividuals <- survivedSims %>%
#   filter(GrowthRate == 4) %>% 
#   group_by(Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()

# lowGrowthTime <- deadSims %>%
#   filter(GrowthRate == 3) %>% 
#   group_by(Time) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# midGrowthTime <- survivedSims %>%
#   filter(GrowthRate == 3.5) %>% 
#   group_by(Time) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# highGrowthTime <- survivedSims %>%
#   filter(GrowthRate == 4) %>% 
#   group_by(Time) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()

# lowGrowthTimeAndInds <- survivedSims %>%
#   filter(GrowthRate == 3) %>% 
#   filter(percentOfFitness != 0) %>% 
#   group_by(Time, Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# midGrowthTimeAndInds <- survivedSims %>%
#   filter(GrowthRate == 3.5) %>% 
#   group_by(Time, Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()
# highGrowthTimeAndInds <- deadSims %>%
#   filter(GrowthRate == 4) %>% 
#   group_by(Time, Individuals) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Time, Individuals, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()

# plot(meanPercentOfFitness~Time,data=lowGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
# plot(meanPercentOfFitness~Time,data=midGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
# plot(meanPercentOfFitness~Time,data=highGrowthTime,ylim=c(0.5,1),xlim=c(0,30))
# plot(meanPercentOfFitness~Individuals,data=lowGrowthIndividuals)
# plot(meanPercentOfFitness~Individuals,data=midGrowthIndividuals)
# plot(meanPercentOfFitness~Individuals,data=highGrowthIndividuals)

# library(lattice)
# wireframe(meanPercentOfFitness ~ Individuals * Time, data=lowGrowthTimeAndInds)
# wireframe(meanPercentOfFitness ~ Individuals * Time, data=midGrowthTimeAndInds)
# wireframe(meanPercentOfFitness ~ Individuals * Time, data=highGrowthTimeAndInds)

# library(ggplot2)
# library(scales)

# colors <- c("#ff0000", 
#             "#ff0000",
#             "#ff0000", 
#             "#ff0000",  
#             "#ff0000", 
#             "#00ffff", 
#             "#ffff00", 
#             "#0000ff", 
#             "#00ff00", 
#             "#000000")

# ggplot(lowGrowthTimeAndInds,aes(x=Individuals,y=Time),) +
#   ggtitle("survived") +
#   geom_tile(aes(fill=meanPercentOfFitness)) +
#   scale_x_continuous(breaks = seq(0, max(lowGrowthTimeAndInds$Individuals), 25),expand = c(0,0)) +
#   scale_y_continuous(breaks = seq(0, 150, 5) ,expand = c(0,0)) +
#   scale_fill_gradientn(name = "", colours = colors ,guide="none") +
#   theme_classic() +
#   labs(fill = "Generations") +
#   theme(legend.title.align=0.5)



# # 
# library(dplyr)
# firstGeneration <- read.csv("./data/out_FitnessInFirstGen_LIFP_1744373846.csv")

# carryingCapacity <- 1000;
# firstGeneration$percentOfFitness <- firstGeneration$Fitness / ( carryingCapacity * firstGeneration$GrowthRate / (carryingCapacity + firstGeneration$GrowthRate * firstGeneration$Females - firstGeneration$Individuals))

# lowGrowthfirstGenerationIF <- firstGeneration %>%
#   filter(GrowthRate == 3) %>% 
#   group_by(Individuals, Females) %>% 
#   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
#   group_by(Individuals, Females, GrowthRate, meanPercentOfFitness) %>% 
#   summarise()

# # lowGrowthfirstGenerationI <- firstGeneration %>%
# #   filter(GrowthRate == 3) %>% 
# #   group_by(Individuals) %>% 
# #   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
# #   group_by(Individuals, GrowthRate, meanPercentOfFitness) %>% 
# #   summarise()

# # lowGrowthfirstGenerationF <- firstGeneration %>%
# #   filter(GrowthRate == 3) %>% 
# #   group_by(Females) %>% 
# #   mutate(meanPercentOfFitness = mean(percentOfFitness))  %>%
# #   group_by(Females, GrowthRate, meanPercentOfFitness) %>% 
# #   summarise()


# library(ggplot2)
# library(scales)

# breaks <- c(0,0.5,0.8,0.9,0.95,1.0)
# colors <- c("#ff0000", "#ff00ff", "#00ff00", "#00ffff", "#0000ff", "#000000")

# # Create a named vector of colors for the scale
# color_scale <- setNames(colors, formatC(breaks, format = "g"))


# ggplot(lowGrowthfirstGenerationIF,aes(x=Individuals,y=Females),) +
#   geom_tile(aes(fill=meanPercentOfFitness)) +
#   scale_x_continuous(breaks = seq(0, max(firstGeneration$Individuals), 5),expand = c(0,0)) +
#   scale_y_continuous(breaks = seq(0, max(firstGeneration$Females), 5) ,expand = c(0,0)) +
#   scale_fill_gradientn(name = "", colours = color_scale, breaks = breaks ,guide="none") +
#   theme_classic() +
#   labs(fill = "Generations") +
#   theme(legend.title.align=0.5)


# plot(meanPercentOfFitness~Females,data=lowGrowthfirstGenerationF,ylim=c(0.5,1),xlim=c(0,50))


# plot(meanPercentOfFitness~Individuals,data=lowGrowthfirstGenerationI,ylim=c(0.5,1),xlim=c(0,65))






library(ggplot2)
library(dplyr)
library(svglite)
fitnesses <- read.csv("./data/out_AlleeEffectFitnessGivenDifferentLoads_LIFPNG_1752936770.csv")
params <- read.csv("./data/out_AlleeEffectFitnessGivenDifferentLoads_LIFPNG_1752936770.csv_params.csv")

colors <- c("#68c24f","#859225", "#81651d","#811d1d","#591414","#a000f0", "#000000")
points <- c(1,3,5,10,15,25,40)
fitnesses$Sterile = params$Sterile[match(fitnesses$Seed,params$Seed)]
fitnesses$PostCompetitionMutationTiming = params$PostCompetitionMutationTiming[match(fitnesses$Seed,params$Seed)]
length(params$Seed)==length(unique(params$Seed))

preSorted <- fitnesses %>% 
    group_by(Seed) %>% 
    arrange(Time)


final_data <- NULL
for (j in 1:length(points)) {
  i <- points[j]
  nthGeneration <- filter(preSorted,Time<=i+1)

  survivedSims <- nthGeneration %>% 
    filter(n()==i) %>%
    summarise(
      finalPercentOfFitness = last(Fitness),
      startSize = first(Individuals),
      startFemales = first(Females),
      secondLastFemales = nth(Females,i-1),
      secondLastSize = nth(Individuals,i-1),
      finalFemales = nth(Females,i),
      finalSize = nth(Individuals,i),
      Sterile,
      PostCompetitionMutationTiming, 
      GrowthRate
    ) %>%
    filter(finalSize!=finalFemales) %>%
    filter(0!=finalFemales)

  final_data_i <- survivedSims %>%
    group_by(startSize, Sterile, PostCompetitionMutationTiming) %>% 
    mutate(meanFinalPercentOfFitness = mean(finalPercentOfFitness))  %>%
    mutate(strongProportion = mean(finalSize>(1/(((GrowthRate-2)/GrowthRate/1000)+1/finalPercentOfFitness/GrowthRate/finalFemales))))  %>%
    mutate(observedGrowthRate = mean((1/(((GrowthRate-2)/GrowthRate/1000)+1/finalPercentOfFitness/GrowthRate/finalFemales))/finalSize))  %>%
    group_by(startSize, meanFinalPercentOfFitness, strongProportion, observedGrowthRate, Sterile, PostCompetitionMutationTiming) %>% 
    summarise() %>%
    mutate(generation = i)
  final_data <- rbind(final_data,final_data_i)
}
final_data$generation <- as.factor(final_data$generation)


svglite("figures/figure_3.svg", width = 8, height = 6)
default_load_data <-filter(final_data, Sterile==0, PostCompetitionMutationTiming==0)
for (j in 1:length(points)) {
  gen <- points[j]
  plot(meanFinalPercentOfFitness~startSize,
    data=filter(default_load_data,generation==gen),
    xlim=c(0,80),
    ylim=c(0.0,1.0),
    col=colors[j],
    pch=16,
    axes=j==1,
    type = "o",
    ylab=if (j==1) "Mean population fitness" else "",
    xlab=if (j==1) "Founding population size (N(0))" else "")
  par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", inset = c(0.05, 0.05),  
  title="Generations",
  legend = points,
  pch=16,
  col=colors)
# abline(h = 2/3, lty = 2,col="red")
dev.off()
# TODO:remove test
svglite("figures/figure_3_test.svg", width = 8, height = 6)
default_load_data <-filter(final_data, Sterile==0, PostCompetitionMutationTiming==0)
for (j in 1:length(points)) {
  gen <- points[j]
  plot(strongProportion~startSize,
    data=filter(default_load_data,generation==gen),
    xlim=c(0,80),
    ylim=c(0.0,1.0),
    col=colors[j],
    pch=16,
    axes=j==1,
    type = "o",
    ylab=if (j==1) "Mean population fitness" else "",
    xlab=if (j==1) "Founding population size (N(0))" else "")
  par(new=TRUE)
}
par(new=FALSE)
legend("topright", inset = c(0.05, 0.05),  
  title="Generations",
  legend = points,
  pch=16,
  col=colors)
dev.off()

# TODO:remove test
svglite("figures/figure_3_test2.svg", width = 8, height = 6)
default_load_data <-filter(final_data, Sterile==0, PostCompetitionMutationTiming==0)
for (j in 1:length(points)) {
  gen <- points[j]
  plot(observedGrowthRate~startSize,
    data=filter(default_load_data,generation==gen),
    xlim=c(0,80),
    ylim=c(0.0,2.0),
    col=colors[j],
    pch=16,
    axes=j==1,
    type = "o",
    ylab=if (j==1) "Mean population fitness" else "",
    xlab=if (j==1) "Founding population size (N(0))" else "")
  par(new=TRUE)
}
par(new=FALSE)
legend("topright", inset = c(0.05, 0.05),  
  title="Generations",
  legend = points,
  pch=16,
  col=colors)
abline(h = 1, lty = 2,col="red")
dev.off()

# svglite("figures/figure_3.svg", width = 8, height = 6)
# points <- c(10,20,40,80)
# default_load_data <-filter(final_data, Sterile==0, PostCompetitionMutationTiming==0)
# for (j in 1:length(points)) {
#   fs <- points[j]
#   plot(meanFinalPercentOfFitness~generation,
#     data=filter(default_load_data,startSize==fs),
#     xlim=c(0,80),
#     ylim=c(0.0,1.0),
#     col=colors[j],
#     pch=16,
#     axes=j==1,
#     type = "o",
#     ylab=if (j==1) "Mean relative fitness" else "",
#     xlab=if (j==1) "Founding population size (N(0))" else "")
#   par(new=TRUE)
# }
# par(new=FALSE)
# legend("bottomright", inset = c(0.05, 0.05),  
#   title="Generations",
#   legend = points,
#   pch=16,
#   col=colors)
# dev.off()

custom_labels <- labeller(
  Sterile = c("0" = "Lethal mutation effect", "1" = "Sterile mutation effect"),
  PostCompetitionMutationTiming = c("0" = "Effect before density-dependent mortality", "1" = "Effect after density-dependent mortality")
)

svglite("figures/figure_S3.svg", width = 8, height = 6)
plot<-ggplot(final_data,aes(x=startSize,y=meanFinalPercentOfFitness,group = generation,colour=generation),) +
  # ggtitle(sprintf("%s gr=%s i=%s",title,gr,i)) +
  geom_line() + geom_point(size = 0.25) +
  scale_color_manual(values = colors) +
  xlab("Founding population size (N(0))") +
  ylab("Mean population fitness") +
  ylim(0.0,1.0) +
  xlim(0,80) +
  theme_light() +
  theme(legend.title.align=0.5) +
  # geom_hline(yintercept = 2/3, linetype = "dashed", color = "red") +
  facet_grid(vars(Sterile), vars(PostCompetitionMutationTiming),labeller = custom_labels) 
print(plot)
dev.off()


# Analysis not included in paper, looks at the extinction rate at each founder size in those populations
params <- read.csv("./data/out_FitnessGivenDifferentLoads30sterile_LIFPNG_1749828832.csv_params.csv")
groupedGraphData <- params %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, PostCompetitionMutationTiming, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, PostCompetitionMutationTiming, RecombinationRate, extinctionRate, count) %>% 
  summarise()
focusedSims <- filter(groupedGraphData,Sterile==0,PostCompetitionMutationTiming==0)
focusedSims2 <- filter(groupedGraphData,Sterile==1,PostCompetitionMutationTiming==0)
plot(extinctionRate~Individuals,data=focusedSims,ylim=c(0,1),ylab="Extinction rate",xlab="Founding population size (N(0))")
par(new=TRUE)
plot(extinctionRate~Individuals,data=focusedSims2,ylim=c(0,1),ylab="Extinction rate",xlab="Founding population size (N(0))",col="#ff0000")



