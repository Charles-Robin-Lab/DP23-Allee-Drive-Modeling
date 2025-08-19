# install.packages("devtools")
# library(devtools)
# install_github("twolodzko/extraDistr")
library(extraDistr)
PMFHomozygotes <- function(N,NMutallele,NHomo) {
    NWTallele <- 2*N-NMutallele
    Nhetero <- NMutallele - 2*NHomo
    NWT<-N-NMutallele+NHomo
    # From Rohlfs & Weir 2008
    exp(lfactorial(N)-lfactorial(NHomo)-lfactorial(Nhetero)-lfactorial(NWT) + lfactorial(NWTallele) + lfactorial(NMutallele)-lfactorial(2*N)+Nhetero*log(2))
}
rHomozygotes <- function(N,Nmuts) {
    # https://stackoverflow.com/questions/42941091/sample-from-custom-distribution-in-r
    r <- runif(1,0,1)
    cdf <- 0
    i <- -1
    while(cdf < r) {
        i <- i+1
        p <- PMFHomozygotes(N,Nmuts,i)
        if (!is.na(p)) {
            cdf <- cdf + p
        }
        if (i>N*2) {
            # I will give $$ if this ever is reached
            # this should only happen if there is a very high r and 
            # floating point error in cdf causing the total AUC to be <1.0 (and by rarest chance <r)
            print("wtf")
            i<-N*2
            break
        }
    }
    i
}

PMFGeneDrift <- function(freq,N,mutCount) {
    # Wright-Fisher Model
    # TODO: I have assumed that changing popsize doesn't affect this formula
    exp(lfactorial(2*N)-lfactorial(mutCount)-lfactorial(2*N-mutCount))*freq^mutCount*(1-freq)^(2*N-mutCount)
}
rGeneDrift <- function(Nmuts,prevN,nextN) {
    r <- runif(1,0,1)
    cdf <- 0
    i <- -1
    while(cdf < r) {
        i <- i+1
        p <- PMFGeneDrift(Nmuts/(2*prevN),nextN,i)
        if (!is.na(p)) {
            cdf <- cdf + p
        }
        if (i>nextN*2) {
            # I will give $$ if this ever is reached
            # this should only happen if there is a very high r and 
            # floating point error in cdf causing the total AUC to be <1.0 (and by rarest chance <r)
            i<-nextN*2
            break
        }
    }
    i
}



# a function for generating samples of S, with W_0 = 1
basicSim <- function(loci, p, Ne, nSamples = 1000) {
    phat <- rbinom(loci * nSamples, size = 2 * Ne, prob = p) # exact binomial draw
    phat <- matrix(phat, nrow = loci) # loci x nSamples matrix
    #   Nhomozygous <- round(Ne*(phat/(2*Ne)^2))
    apply(phat, 2, function(x) {
        Hi <- (x / (2 * Ne))^2
        prod(1 - Hi)
    })
}
round_stochastic <- function(x) {
    ifelse(runif(length(x),0,1) > x %% 1,floor(x),ceiling(x))
}

genotypeSim <- function(loci, p, Ne, birthRate, carryingCapacity, mortalityRate, nSamples = 1000) {
    out = rep(Ne, nSamples)
    for (i in 1:nSamples) {
        currentN <- round_stochastic(birthRate*Ne)
        currentMutCounts <- round_stochastic(birthRate*rbinom(loci, size = 2 * Ne, prob = p))
        # currentMutCounts<-rGeneDrift(currentMutCounts[j],Ne,nextN)
        
        # initialise with correct size 
        currentGenotypeCounts <- matrix(0, 3, loci)
        # initialise genotype dists
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
        }
        # density-dependent mortality
        # currentN <- 1 / (mortalityRate+1/currentN)
        out[i] <- currentN
    }
    out
}

timeTilNextDeath <- function(rate,roll) {
    -log(1-roll)/rate
}
rMortalityN <- function(u,N) {
    r <- runif(1,0,1)
    passed_time <- timeTilNextDeath(u*N^2,r)
    while (passed_time < 1 && N>0) {
        # print(passed_time)
        # print(N)
        N <- N-1
        r <- runif(1,0,1)
        passed_time <- passed_time + timeTilNextDeath(u*N^2,r)
    }
    N
}


recursiveSim <- function(loci, p, Ne, birthRate, carryingCapacity, mortalityRate, nSamples = 1000) {
    out = rep(Ne, nSamples)
    for (i in 1:nSamples) {
        #  initial params
        currentN <- Ne
        currentMutCounts <- rbinom(loci, size = 2 * Ne, prob = p)
        # initialise with correct size 
        currentGenotypeCounts <- matrix(0, 3, loci)
        generations <- 0
        while(currentN > 0 && currentN < carryingCapacity && generations <= 2000) {
            generations <- generations + 1
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
            # density-dependent mortality
            currentN <- rMortalityN(mortalityRate,currentN)
            # update all loci genotype counts after all the eliminations
            for(j in 1:loci) {
                currentGenotypeCounts[,j] <- rmvhyper(1,currentGenotypeCounts[,j],currentN)
                currentMutCounts[j] <- currentGenotypeCounts[2,j]
            }
        }
        out[i] <- if (currentN<=0) "EXTINCT" else  if (generations>=2000) "LOADED_SURVIVAL" else "SURVIVAL"
    }
    out
}