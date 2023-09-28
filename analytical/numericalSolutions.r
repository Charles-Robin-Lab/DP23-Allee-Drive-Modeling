library(extraDistr)
# library(devtools)
# install_github("twolodzko/extraDistr")
PMFHomozygotes <- function(N,NMutallele,NHomo) {
    NWTallele <- 2*N-NMutallele
    Nhetero <- NMutallele - 2*NHomo
    NWT<-N-NMutallele+NHomo
    # From Rohlfs & Weir 2008
    exp(lfactorial(N)-lfactorial(NHomo)-lfactorial(Nhetero)-lfactorial(NWT) + lfactorial(NWTallele)-lfactorial(NMutallele)-lfactorial(2*N)) * 2^Nhetero
}
rHomozygotes <- function(N,Nmuts) {
    # https://stackoverflow.com/questions/42941091/sample-from-custom-distribution-in-r
    r <- runif(1,0,1)
    cdf <- 0
    i <- -1
    while(cdf < r) {
        i <- i+1
        p <- PMFHomozygotes(N,Nmuts,i)
        
        # if (length(p)==0) {
        #     print("Homo")
        #     print(Nmuts)
        #     print(N)
        # }
        if (!is.na(p)) {
           cdf <- cdf + p
        }
        if (i>N*2) {
            # print("Error! too much growth")
            i<-round(N*(Nmuts/(2*N))^2)
            break
        }
    }
    i
}

PMFGeneDrift <- function(freq,N,mutCount) {
    # Wright-Fisher Model
    # I have assumed that changing popsize doesn't affect this formula
    exp(lfactorial(2*N)-lfactorial(mutCount)-lfactorial(2*N-mutCount))*freq^mutCount*(1-freq)^(2*N-mutCount)
}
rGeneDrift <- function(Nmuts,prevN,nextN) {
    r <- runif(1,0,1)
    cdf <- 0
    i <- -1
    while(cdf < r) {
        i <- i+1
        p <- PMFGeneDrift(Nmuts/(2*prevN),nextN,i)
        # if (length(p)==0) {
        #     print("GD")
            # print(cdf)
            # print(p)
        # }
        if (!is.na(p)) {
           cdf <- cdf + p
        }
        if (i>nextN*2) {
            # print("Error! too much growth")
            i<-Nmuts
            break
        }
    }
    i
}



# a function for generating samples of W, with W_0 = 1
W <- function(loci, p, Ne, nSamples = 1000) {
    phat <- rbinom(loci * nSamples, size = 2 * Ne, prob = p) # exact binomial draw
    phat <- matrix(phat, nrow = loci) # n x nSamples matrix
    #   sample binomial amount of homozygotes
    #   Nhomozygous <- round(Ne*(phat/(2*Ne)^2))
    apply(phat, 2, function(x) {
        Nhomozygous <- Ne * (x / (2 * Ne))^2
        prod(1 - Nhomozygous / Ne)
    })
}

extendedSim <- function(n, p, Ne, nSamples = 1000) {
    out = rep(Ne, nSamples)
    for (i in 1:nSamples) {
        currentN <- Ne
       for (j in 1:loci) {
            NMutallele <- rbinom(1, size = 2 * currentN, prob = p)
            Nhomozygous <- rHomozygotes(currentN,NMutallele)
            currentN <- currentN - Nhomozygous
       }
        out[i] <- currentN
    }
    out
}

extendedSimRecurse <- function(loci, p, Ne, birthRate, carryingCapacity, nSamples = 1000) {
    out = rep(Ne, nSamples)
    for (i in 1:nSamples) {
        #  initial params
        currentN <- Ne
        currentMutCounts <- rbinom(loci, size = 2 * Ne, prob = p)
        # initialise with correct size 
        currentGenotypeCounts <- matrix(0, 3, loci)
       while(currentN > 0 && currentN < 80/birthRate) {
            # Genotype dists
            for (j in 1:loci) {
               currentGenotypeCounts[1,j] <- rHomozygotes(currentN,currentMutCounts[j])
               currentGenotypeCounts[2,j] <- currentMutCounts[j] - 2*currentGenotypeCounts[1,j]
               currentGenotypeCounts[3,j] <- currentN - currentGenotypeCounts[2,j] - currentGenotypeCounts[1,j]
            }
            # eliminate Homozygotes
            for (j in 1:loci) {
                    # update genotype counts based on the already eliminated individuals
                    currentGenotypeCounts[,j] <- rmvhyper(1,currentGenotypeCounts[,j],currentN)
                    currentN <- currentN - currentGenotypeCounts[1,j]
                    currentGenotypeCounts[1,j] <- 0
            }
            # update all loci genotype counts after all the eliminations
            for(j in 1:loci) {
                currentGenotypeCounts[,j] <- rmvhyper(1,currentGenotypeCounts[,j],currentN)
                currentMutCounts[j] <- currentGenotypeCounts[2,j]
            }
            # reproduce and drift
            nextN <- rpois(1,birthRate*currentN)
            for (j in 1:loci) {
                currentMutCounts[j] <- rGeneDrift(currentMutCounts[j],currentN,nextN)
            }
            currentN <- nextN
       }
        out[i] <- if (currentN<=0) "EXTINCT" else "SURVIVAL"
    }
    out
}
extendedSimRecurse2 <- function(loci, p, Ne, birthRate, carryingCapacity, nSamples = 1000) {
    out = rep(Ne, nSamples)
    for (i in 1:nSamples) {
        #  initial params
        currentN <- Ne
        currentMutCounts <- rbinom(loci, size = 2 * Ne, prob = p)
        # initialise with correct size 
        currentGenotypeCounts <- matrix(0, 3, loci)
       while(currentN > 0 && currentN < 80/birthRate) {
            # reproduce and drift
            nextN <- rpois(1,birthRate*currentN)
            for (j in 1:loci) {
                currentMutCounts[j] <- rGeneDrift(currentMutCounts[j],currentN,nextN)
            }
            currentN <- nextN
            # Genotype dists
            for (j in 1:loci) {
               currentGenotypeCounts[1,j] <- rHomozygotes(currentN,currentMutCounts[j])
               currentGenotypeCounts[2,j] <- currentMutCounts[j] - 2*currentGenotypeCounts[1,j]
               currentGenotypeCounts[3,j] <- currentN - currentGenotypeCounts[2,j] - currentGenotypeCounts[1,j]
            }
            # eliminate Homozygotes
            for (j in 1:loci) {
                    # update genotype counts based on the already eliminated individuals
                    currentGenotypeCounts[,j] <- rmvhyper(1,currentGenotypeCounts[,j],currentN)
                    currentN <- currentN - currentGenotypeCounts[1,j]
                    currentGenotypeCounts[1,j] <- 0
            }
            # update all loci genotype counts after all the eliminations
            for(j in 1:loci) {
                currentGenotypeCounts[,j] <- rmvhyper(1,currentGenotypeCounts[,j],currentN)
                currentMutCounts[j] <- currentGenotypeCounts[2,j]
            }
       }
        out[i] <- if (currentN<=0) "EXTINCT" else "SURVIVAL"
    }
    out
}




# some example parameters...
loci <- 100 # number of loci
freqRange <- seq(0.0, 0.3, by = 0.0025) # mean frequency in ancestral population
Ne <- 25
birthRate <- 3.0
carryingCapacity <- 1000
survivalW <- list()
for (averageMutFreq in freqRange) {
    mutFreqs <- rep(averageMutFreq,loci)
    fitnessDecreases <- W(loci, mutFreqs, Ne)
    survivalW <- append(survivalW, length(fitnessDecreases[fitnessDecreases > 1 / birthRate]) / length(fitnessDecreases))
}
par(new=TRUE)
plot(freqRange, survivalW,col="red",ylim=c(0,1),pch=2,xlab="",ylab="")

survivalSingle <- list()
for (averageMutFreq in freqRange) {
    fitnessDecreases <- extendedSim(loci, averageMutFreq, Ne)/Ne
    survivalSingle <- append(survivalSingle, length(fitnessDecreases[fitnessDecreases > 1 / birthRate]) / length(fitnessDecreases))
}
par(new=TRUE)
plot(freqRange, survival,col="green",ylim=c(0,1),pch=2,xlab="",ylab="")

survivalRecurse <- list()
for (averageMutFreq in freqRange) {
    outcomeGroup <- extendedSimRecurse2(loci, averageMutFreq, Ne, birthRate, carryingCapacity)
    # fitnessDecreases <- W(loci, mutFreqs, Ne)
    survivalRecurse <- append(survivalRecurse, length(outcomeGroup[outcomeGroup != "EXTINCT"]) / length(outcomeGroup))
}

par(new=TRUE)
plot(head(freqRange,length(survivalRecurse)), survivalRecurse,col="blue",ylim=c(0,1),xlim=c(0,0.3),pch=0,xlab="",ylab="")


# Pulling slim apart:
# a- initial sampling bias
# s- binomial sexing
# a- lethality
# a- variance in reproductive output (poisson of lambda)
# a- genetic drift (allele frequency)
# a- HW drift (genotype frequency)
# x- Distribution of alleles within individuals