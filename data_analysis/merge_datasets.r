# Tool for merging datasets that have been split to achieve further parralelisation.
# aditionally verifies the 400 replicates are present for each paramater combination, there can be extra if compute nodes crash and restart their job
library(dplyr)



dataset0<- "./data/out_bestParameterSpaceFirstQuarter_LIFPNG_1752545199.csv"
dataset1<- "./data/out_bestParameterSpaceSecondQuarter_LIFPNG_1752545204.csv"
dataset2<- "./data/out_bestParameterSpaceFourthQuarter_LIFPNG_1752545194.csv"
dataset3<- "./data/out_bestParameterSpaceThirdQuarter_LIFPNG_1752545209.csv"
outputfile<- "E:/out_bestParameterSpace_LIFPNG_1752545194.csv"

groupedGraphData <- bind_rows(read.csv(dataset0),read.csv(dataset1),read.csv(dataset2),read.csv(dataset3)) %>%
    group_by(
        Individuals,
        GrowthRate,
        RecombinationRate,
        ChromosomeCount,
        MaxGenerations,
        CarryingCapacity,
        MutationFrequency,
        MutationCount,
        Sterile,
        PostCompetitionMutationTiming,
        Xlinked,
        FemaleOnlyEffect)
verify <- groupedGraphData %>%
        mutate (count=n()) %>%
        group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, FemaleOnlyEffect, count) %>% 
        summarise()
if (any(verify$count!=400)) {
    hist(verify$count)
} else {
    print("all data verified")
}

cleaned <- groupedGraphData %>% 
    slice_head(n= 400) %>%
    ungroup()

write.csv(cleaned, outputfile, row.names = FALSE)