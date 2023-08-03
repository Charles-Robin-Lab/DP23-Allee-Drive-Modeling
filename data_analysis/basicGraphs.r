library(dplyr)

groupedData <- read.csv("./data/out.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()
autosomalData <- groupedData[groupedData$Xlinked==0,]
groupedGraphData <- read.csv("./data/graphData.csv") %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked) %>% 
  mutate(count = n()) %>%
  mutate(survivalRate = sum(Result == "SURVIVED") / count) %>%
  group_by(MutationFrequency, MutationCount, Individuals, GrowthRate, Sterile, Xlinked, survivalRate,count) %>% 
  summarise()
autosomalGraphData <- groupedGraphData[groupedGraphData$Xlinked==0,]


dataSlice2d.MutationFrequency <- filter(autosomalData, MutationCount==61, Individuals == 22, GrowthRate == 2, Sterile==1)
dataSlice2d.MutationCount <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==1)
dataSlice2d.Lethal <- filter(autosomalData, MutationFrequency==0.11, Individuals == 22, GrowthRate == 2, Sterile==0)
dataSlice2d.Individuals <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, GrowthRate == 2, Sterile==1)
dataSlice2d.GrowthRate <- filter(autosomalData, MutationFrequency==0.11, MutationCount==61, Individuals == 22, Sterile==1)




# MutationCount==60, MutationFrequency==0.10, Individuals == 20, GrowthRate == 1, Sterile==1

dataSlice2d.MutationFrequency <- filter(autosomalGraphData, MutationCount==60, Individuals == 20, GrowthRate == 1, Sterile==1)
dataSlice2d.MutationCount <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 1, Sterile==1)
dataSlice2d.Lethal <- filter(autosomalGraphData, MutationFrequency==0.10, Individuals == 20, GrowthRate == 2, Sterile==0)
dataSlice2d.Individuals <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, GrowthRate == 1, Sterile==1)
dataSlice2d.GrowthRate <- filter(autosomalGraphData, MutationFrequency==0.10, MutationCount==60, Individuals == 20, Sterile==1)



plot(survivalRate~MutationFrequency,data=dataSlice2d.MutationFrequency,ylim=c(0,1))
plot(survivalRate~MutationCount,data=dataSlice2d.MutationCount,ylim=c(0,1))
plot(survivalRate~Individuals,data=dataSlice2d.Individuals,ylim=c(0,1))
plot(survivalRate~GrowthRate,data=dataSlice2d.GrowthRate,ylim=c(0,1))




dataSlice3d.countfreq <- filter(groupedData, Individuals == 32, GrowthRate == 2, Sterile==1, Xlinked==1)

dataSlice3d.countInd <- filter(groupedData, MutationFrequency == 0.11, GrowthRate == 2, Sterile==1, Xlinked==1)

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
