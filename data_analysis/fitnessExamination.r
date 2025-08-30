# # entire sim fitness trends
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


# Analysis not included in paper, looks at the extinction proportion at each founder size in the above populations
params <- read.csv("./data/out_AlleeEffectFitnessGivenDifferentLoads_LIFPNG_1752936770.csv_params.csv")
groupedGraphData <- params %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, PostCompetitionMutationTiming, RecombinationRate) %>% 
  mutate(count = n()) %>%
  mutate(extinctionRate = sum(Result == "EXTINCT") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, PostCompetitionMutationTiming, RecombinationRate, extinctionRate, count) %>% 
  summarise()
focusedSims <- filter(groupedGraphData,Sterile==0,PostCompetitionMutationTiming==0)
focusedSims2 <- filter(groupedGraphData,Sterile==1,PostCompetitionMutationTiming==0)
plot(extinctionRate~Individuals,data=focusedSims,ylim=c(0,1),ylab="Proportion extinct",xlab="Founding population size (N(0))")
par(new=TRUE)
plot(extinctionRate~Individuals,data=focusedSims2,ylim=c(0,1),ylab="Proportion extinct",xlab="Founding population size (N(0))",col="#ff0000")



