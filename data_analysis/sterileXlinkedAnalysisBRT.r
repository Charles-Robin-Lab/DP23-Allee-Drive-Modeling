# A script to examine how extinction rate of sterile vs lethal effect varies through parameter space
  # Exploratory analysis using boosted regression trees
  # response variable is the normalised difference in extinction probability
  # explanatory variables are parameter axes

library(tidyverse)

########### Prepare the data ###########

# load dataset
simulatedDataset <- "./data/out_parameterSpace_LIFP_1715053595.csv"
d <- read.csv(simulatedDataset)

# Summarise to give extinction probability
d <- d %>%
  # remove non-varying parameters
  select(-RecombinationRate, -ChromosomeCount, -MaxGenerations, -CarryingCapacity, -FemaleOnlyEffect, -Seed) %>%
  # group by remaining parameter axes
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  # calculate extinction probability for each part of parameter space
  summarise(count = n(), ExtinctionRate = sum(Result == "EXTINCT") / count)

# Calculate difference between Sterile/lethal extinction for each part of parameter space
modelData <- d %>%
  # pivot wider
  pivot_wider(names_from = Sterile, values_from = ExtinctionRate, names_prefix = "Sterile") %>%
  # calculate difference in extinction probability on the logit scale
  mutate(logitLethal = log(Sterile0/(1-Sterile0)), logitSterile = log(Sterile1/(1-Sterile1)), diffExtinctionRate = logitSterile-logitLethal) %>%
  # remove situations where extinction probability was 1 or 0
  filter(is.finite(diffExtinctionRate)) %>%
  ungroup() %>%
  as.data.frame() # else dismo cries
  

boxplot(modelData$diffExtinctionRate)

####### Use diffExtinctionRate as response variable in BRT #########

# load dismo here to avoid conflict with select() in dplyr
library(dismo)

# fit a model
mod <- gbm.step(data = modelData, gbm.x = 1:5, gbm.y = 11, 
                family = "gaussian", 
                tree.complexity = 4,
                learning.rate = 0.01, 
                bag.fraction = 0.5)

# report on variable importance
summary(mod)

# make a plot of marginal predictor functions (omitting x-linked which doesn't do much)
pdf("out/BRT-marginals.pdf")
  gbm.plot(mod, n.plots = 4, 
           plot.layout = c(2,2), 
           rug = FALSE, 
           y.label = "Difference in extinction rate",
           write.title = FALSE)
dev.off()

# Examine interactions
modInt <- gbm.interactions(mod)
modInt$rank.list
modInt$interactions

# plot the most interesting two interactions
pdf("out/BRT-interactions.pdf")
  gbm.perspec(mod, 2, 1, z.range = c(-1, 5))
  gbm.perspec(mod, 4, 3, z.range = c(-1, 9))
dev.off()
