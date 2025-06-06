maxGrowthRate=1.0;
carryingCapacity=1000;
N=10
intitialMutFreq=0.10
mutatedGenes=100;

print(N)
mutFreq = rep(intitialMutFreq,mutatedGenes);
# print(mutFreq[1])
while (N>0 && N<0.5*carryingCapacity) {
    # kill off things
    for (i in 1:mutatedGenes) {
        viableN=(1-mutFreq[i]^2)*N
        nMutCount = 2*mutFreq[i]*(1-mutFreq[i])*N
        nMutFreq=nMutCount/(viableN*2)
        N=viableN
        
        mutFreq[i]=nMutFreq
    }
    N=round(N)
    mutFreq=round(mutFreq*2*N)/(2*N)
    # reproduce
    nN = N + maxGrowthRate*N*(1-N/carryingCapacity) ;
    nN=round(nN)
    mutFreq=round(mutFreq*2*nN)/(2*nN)

    # t+1
    N = nN;
    print(N)
    # print(mutFreq[1])
}


# Pulling slim apart:
# - initial sampling bias
# - binomial sexing
# - lethality
# - variance in reproductive output (poisson of lambda)
# - genetic drift (allele frequency)
# - HW drift (genotype frequency)


# t=0
# mutFreq=intitialMutFreq
# print(N)
# while (N<0.5*carryingCapacity) {
#     N=(1-(1/(t+1/mutFreq))^2)^mutatedGenes*N
#     N= N + maxGrowthRate*N*(1-N/carryingCapacity) ;
#     t=t+1
#     print(N)
# }


# discretising:
# - are these males or females

# error bars
# variability in: 
# - sex ratio
# - allele frequency
        # genedriftsd = sqrt(nMutFreq*(1-nMutFreq)/(2*N))
# - ?HW equilibrium?
