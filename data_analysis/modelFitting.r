library(dplyr)
library(svglite)
library("gridExtra")
library(ggplot2)

options(width=2000)




dataset <-"./data/out_LoadEquivalenceParameterSpace_LIFPNG_1753296009.csv"


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

svglite("figures/figure_4.svg", width = 13*dotpointscale, height = 18*dotpointscale)
par(mar=c(5,4,2,2)+0.1)
ggplot(groupedGraphData,aes(x=MutationFrequency,y=MutationCount),) +
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
    print(sprintf("b=%d",gr))
    print(sprintf("N(0)=%d",i))
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
    
    print("-----------")
    print(sprintf("b=%d",gr))
    print(sprintf("N(0)=%d",i))
    step(model,direction='both')
    print(sprintf("b=%d",gr))
    print(sprintf("N(0)=%d",i))
  }
}