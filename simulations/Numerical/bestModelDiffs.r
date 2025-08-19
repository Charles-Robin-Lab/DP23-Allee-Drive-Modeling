source("./slimulations/Numerical/numericalModels.r")
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (JobID).n", call.=FALSE)
}
jobid <- args[1]
lociRange <- seq(10, 250, by = 40)
freqRange <- seq(0.0, 0.15, by = 0.01)
N0Range <- seq(10, 100, by = 10)
birthRateRange <- seq(1.0, 3.0, by = 0.25)
carryingCapacity <- 1000

write(paste("loci,freq,N0,birthRate,carryingCapacity",sep = ""),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""))
write(paste(paste(lociRange,collapse=' '),paste(freqRange,collapse=' '),paste(N0Range,collapse=' '),paste(birthRateRange,collapse=' '),carryingCapacity,sep = ","),paste("./data/",as.character(jobid),"_job_params.csv",sep = ""),append=TRUE)


for (iloci in lociRange) {
    for (iaverageMutFreq in freqRange) {
        for (iN0 in N0Range) {
            survivalRecurse <- list()
            for (ibirthRate in birthRateRange) {
                imortalityRate <- (ibirthRate - 1)/(ibirthRate*carryingCapacity)
                # Multigeneration
                outcomeGroup <- recursiveSim(iloci, iaverageMutFreq, iN0, ibirthRate, carryingCapacity, imortalityRate, 1)
                survivalRecurse <- append(survivalRecurse, outcomeGroup)
            }
            write.table(data.frame(rep(iloci,length(survivalRecurse)),rep(iaverageMutFreq,length(survivalRecurse)),rep(iN0,length(survivalRecurse)),birthRateRange,unlist(survivalRecurse)), 
                        file = paste("./data/",as.character(jobid),"_model3.csv",sep = ""), sep = ",", 
                        append = TRUE, quote = FALSE, 
                        col.names = FALSE, row.names = FALSE)
        }
    }
}


