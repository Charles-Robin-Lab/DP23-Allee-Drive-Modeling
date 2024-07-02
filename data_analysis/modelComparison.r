
library(svglite)
library(ggplot2)
library(dplyr)
freqRange <- read.csv("./data/analytical_model1.csv")$freqRange
survivalW <- read.csv("./data/analytical_model1.csv")$unlist.survivalW.
survivalSingle <- read.csv("./data/analytical_model2_2.csv")$unlist.survivalSingle.
# survivalRecurse <- read.csv("./data/analytical_model3.csv")$unlist.survivalRecurse.
survivalRecurse <- read.csv("./data/analytical_model3_rerun.csv")$unlist.survivalRecurse.
dataSlice2d.MutationFrequency <- read.csv("./data/modelslim.csv")
A = data.frame(x = freqRange, y = survivalW)
B = data.frame(x = freqRange, y = survivalSingle)
C = data.frame(x = freqRange, y = survivalRecurse)
D = data.frame(x = dataSlice2d.MutationFrequency$MutationFrequency, y = dataSlice2d.MutationFrequency$extinctionProbability)

A = data.frame(x = rnorm(10),y=rnorm(10))
B = data.frame(x = rnorm(10),y=rnorm(10))
ggplot(A,aes(x,y)) +geom_point() +geom_point(data=B,colour='red') + xlim(0, 10)


svglite("figures/figure_1.svg", width = 8, height = 6)
plot(freqRange, survivalW,col="red",ylim=c(0,1),xlim=c(0.0,0.2),pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(freqRange, survivalSingle,col="green",ylim=c(0,1),xlim=c(0.0,0.2),pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(freqRange,survivalRecurse,col="blue",ylim=c(0,1),xlim=c(0.0,0.2),pch=2,xlab="",ylab="",axes = FALSE)
par(new=TRUE)
plot(extinctionProbability~MutationFrequency,data=filter(dataSlice2d.MutationFrequency),col="black",ylim=c(0,1),xlim=c(0.0,0.2),xlab="Deleterious recessive frequency",ylab="Extinction probability")
dev.off()