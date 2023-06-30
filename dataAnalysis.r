data <- read.csv("./out.csv")
# summary(lm(Result ~ MutationFrequency*MutationCount,data=data))



library(dplyr)

# Assuming your data frame is named 'data'
# Create a new column 'count' that counts the occurrences of each combination of a, b, and c
groupedData <- data %>%
  dplyr::group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>%
  dplyr::mutate(count = n())

# # Create a new column 'T_percentage' that calculates the percentage of T results for each combination of a, b, and c
groupedData <- groupedData %>%
  dplyr::group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>%
  dplyr::mutate(T_percentage = sum(Result == "SURVIVED") / count * 100)


data1 <- filter(groupedData, MutationFrequency==0.16, Individuals == 32, GrowthRate == 2, Sterile==1, Xlinked==0)
plot(data1$MutationCount,data1$T_percentage)

# sapply(data1$MutationFrequency,head,100)

# View the modified data frame
print(data1)
# print(filter(data, Time>5000, Result == "EXTINCT"))


