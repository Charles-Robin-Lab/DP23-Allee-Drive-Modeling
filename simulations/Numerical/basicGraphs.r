source(numericalModels.r)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (JobID).n", call.=FALSE)
}
jobid <- args[1]
loci <- 200 # number of loci #100
lociRange <- seq(2, 250, by = 2)
freq <- 0.04
averageMutFreq <- freq
freqRange <- seq(0.0, 0.1, by = 0.00125) # mean frequency in ancestral population
N0 <- 25 # 10
N0Range <- seq(2, 120, by = 1)
birthRate <- 1.5 #5
birthRateRange <- seq(1.0, 2.0, by = 0.02)
carryingCapacity <- 1000
mortalityRate <- (birthRate - 1)/(birthRate*carryingCapacity)

write(paste("loci,freq,N0,birthRate,carryingCapacity",sep = ""),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""))
write(paste(paste(lociRange,collapse=' '),freq,N0,birthRate,carryingCapacity,sep = ","),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""),append=TRUE)
write(paste(loci,paste(freqRange,collapse=' '),N0,birthRate,carryingCapacity,sep = ","),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""),append=TRUE)
write(paste(loci,freq,paste(N0Range,collapse=' '),birthRate,carryingCapacity,sep = ","),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""),append=TRUE)
write(paste(loci,freq,N0,paste(birthRateRange,collapse=' '),carryingCapacity,sep = ","),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""),append=TRUE)


survivalW <- list()
survivalSingle <- list()
survivalRecurse <- list()
for (iloci in lociRange) {
    # Basic
    mutFreqs <- rep(averageMutFreq,iloci)
    lambda <- (1/(mortalityRate + 1/(birthRate*N0*basicSim(iloci, mutFreqs, N0, 1))))/N0
    survivalW <- append(survivalW, if (lambda < 1) "EXTINCT" else "SURVIVAL")    
    # Random Genotype
    lambda <- (genotypeSim(iloci, averageMutFreq, N0, birthRate, carryingCapacity, mortalityRate, 1))/N0
    survivalSingle <- append(survivalSingle, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Multigeneration
    outcomeGroup <- recursiveSim(iloci, averageMutFreq, N0, birthRate, carryingCapacity, mortalityRate, 1)
    survivalRecurse <- append(survivalRecurse, outcomeGroup)
}
write.table(data.frame(lociRange,rep(freq,length(survivalW)),rep(N0,length(survivalW)),rep(birthRate,length(survivalW)),unlist(survivalW)), 
            file = paste("./data/",as.character(jobid),"_model1.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(lociRange,rep(freq,length(survivalSingle)),rep(N0,length(survivalSingle)),rep(birthRate,length(survivalSingle)),unlist(survivalSingle)), 
            file = paste("./data/",as.character(jobid),"_model2.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(lociRange,rep(freq,length(survivalRecurse)),rep(N0,length(survivalRecurse)),rep(birthRate,length(survivalRecurse)),unlist(survivalRecurse)), 
            file = paste("./data/",as.character(jobid),"_model3.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

survivalW <- list()
survivalSingle <- list()
survivalRecurse <- list()
for (iaverageMutFreq in freqRange) {
    # Basic
    mutFreqs <- rep(iaverageMutFreq,loci)
    lambda <- (1/(mortalityRate + 1/(birthRate*N0*basicSim(loci, mutFreqs, N0, 1))))/N0
    survivalW <- append(survivalW, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Random Genotype
    lambda <- (genotypeSim(loci, iaverageMutFreq, N0, birthRate, carryingCapacity, mortalityRate, 1))/N0
    survivalSingle <- append(survivalSingle, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Multigeneration
    outcomeGroup <- recursiveSim(loci, iaverageMutFreq, N0, birthRate, carryingCapacity, mortalityRate, 1)
    survivalRecurse <- append(survivalRecurse, outcomeGroup)
}
write.table(data.frame(rep(loci,length(survivalW)),freqRange,rep(N0,length(survivalW)),rep(birthRate,length(survivalW)),unlist(survivalW)), 
            file = paste("./data/",as.character(jobid),"_model1.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalSingle)),freqRange,rep(N0,length(survivalSingle)),rep(birthRate,length(survivalSingle)),unlist(survivalSingle)), 
            file = paste("./data/",as.character(jobid),"_model2.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalRecurse)),freqRange,rep(N0,length(survivalRecurse)),rep(birthRate,length(survivalRecurse)),unlist(survivalRecurse)), 
            file = paste("./data/",as.character(jobid),"_model3.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)


survivalW <- list()
survivalSingle <- list()
survivalRecurse <- list()
for (iN0 in N0Range) {
    # Basic
    mutFreqs <- rep(averageMutFreq,loci)
    lambda <- (1/(mortalityRate + 1/(birthRate*iN0*basicSim(loci, mutFreqs, iN0, 1))))/N0
    survivalW <- append(survivalW, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Random Genotype
    lambda <- (genotypeSim(loci, averageMutFreq, iN0, birthRate, carryingCapacity, mortalityRate, 1))/iN0
    survivalSingle <- append(survivalSingle, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Multigeneration
    outcomeGroup <- recursiveSim(loci, averageMutFreq, iN0, birthRate, carryingCapacity, mortalityRate, 1)
    survivalRecurse <- append(survivalRecurse, outcomeGroup)
}
write.table(data.frame(rep(loci,length(survivalW)),rep(freq,length(survivalW)),N0Range,rep(birthRate,length(survivalW)),unlist(survivalW)), 
            file = paste("./data/",as.character(jobid),"_model1.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalSingle)),rep(freq,length(survivalSingle)),N0Range,rep(birthRate,length(survivalSingle)),unlist(survivalSingle)), 
            file = paste("./data/",as.character(jobid),"_model2.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalRecurse)),rep(freq,length(survivalRecurse)),N0Range,rep(birthRate,length(survivalRecurse)),unlist(survivalRecurse)), 
            file = paste("./data/",as.character(jobid),"_model3.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)


survivalW <- list()
survivalSingle <- list()
survivalRecurse <- list()
for (ibirthRate in birthRateRange) {
    imortalityRate <- (ibirthRate - 1)/(ibirthRate*carryingCapacity)
    # Basic
    mutFreqs <- rep(averageMutFreq,loci)
    lambda <- (1/(imortalityRate + 1/(ibirthRate*N0*basicSim(loci, mutFreqs, N0, 1))))/N0
    survivalW <- append(survivalW, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Random Genotype
    lambda <- (genotypeSim(loci, averageMutFreq, N0, ibirthRate, carryingCapacity, imortalityRate, 1))/N0
    survivalSingle <- append(survivalSingle, if (lambda < 1) "EXTINCT" else "SURVIVAL")
    # Multigeneration  
    outcomeGroup <- recursiveSim(loci, averageMutFreq, N0, ibirthRate, carryingCapacity, imortalityRate, 1)
    survivalRecurse <- append(survivalRecurse, outcomeGroup)  
}
write.table(data.frame(rep(loci,length(survivalW)),rep(freq,length(survivalW)),rep(N0,length(survivalW)),birthRateRange,unlist(survivalW)), 
            file = paste("./data/",as.character(jobid),"_model1.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalSingle)),rep(freq,length(survivalSingle)),rep(N0,length(survivalSingle)),birthRateRange,unlist(survivalSingle)), 
            file = paste("./data/",as.character(jobid),"_model2.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
write.table(data.frame(rep(loci,length(survivalRecurse)),rep(freq,length(survivalRecurse)),rep(N0,length(survivalRecurse)),birthRateRange,unlist(survivalRecurse)), 
            file = paste("./data/",as.character(jobid),"_model3.csv",sep = ""), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)



