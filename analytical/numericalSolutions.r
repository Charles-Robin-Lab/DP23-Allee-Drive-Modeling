
PMFHomozygote <- function(N,NMutallele,NHomo) {
    # NMutallele <- round(2*N*p)
    NWTallele <- 2*N-NMutallele
    Nhetero <- NMutallele - 2*NHomo
    NWT<-N-NMutallele+NHomo
    # if (any(list(N,NMutallele,NHomo,NWT,NWTallele,Nhetero)<0)) {
    #     print(list(N,NMutallele,NHomo,NWT,NWTallele,Nhetero))
    # }
    # From Rohlis & Weir 2008
    (factorial(N)/factorial(NHomo)/factorial(Nhetero)/factorial(NWT)) * 2^Nhetero * factorial(NWTallele)* factorial(NMutallele)/factorial(2*N)
}
rHomozygotes <- function(N,Nmuts) {
    r <- runif(1,0,1)
    cdf <- 0
    i <- -1
    while(cdf < r){
        i <- i+1
        p <- PMFHomozygote(N,Nmuts,i)
        if (!is.na(p)) {
           cdf <- cdf + p
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
    out = rep(N, nSamples)
    for (i in 1:nSamples) {
        currentN <- Ne
       for (j in 1:loci) {
            NMutallele <- rbinom(1, size = 2 * currentN, prob = p)
            Nhomozygous <- rHomozygotes(currentN,NMutallele)
            currentN <- currentN - Nhomozygous
       }
        out[i-1] <- currentN
    }
    out
}





# some example parameters...
loci <- 100 # number of loci
freqRange <- seq(0, 0.3, by = 0.0025) # mean frequency in ancestral population
Ne <- 25
birthRate <- 2.0
survival <- list()
for (averageMutFreq in freqRange) {
    mutFreqs <- rep(averageMutFreq, loci) # vector of p's
    fitnessDecreases <- extendedSim(loci, mutFreqs, Ne)/Ne
    # fitnessDecreases <- W(loci, mutFreqs, Ne)
    survival <- append(survival, length(fitnessDecreases[fitnessDecreases > 1 / birthRate]) / length(fitnessDecreases))
}
plot(freqRange, survival)

# hOut <- hist(W(p, Ne),
#              freq = FALSE,
#              main = "",
#              xlab = "Fitness")
# text(x = hOut$mids[1], y = max(hOut$density), labels = paste("Ne = ", Ne))
