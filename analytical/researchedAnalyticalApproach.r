# based off Yoshinari Tanaka 2000
count = 0;
replicates = 100;
for (i in 1:replicates) {
    averageMutFreq=0.16;
    maxGrowthRate=1.0;
    N=32;
    mutations=81;
    selectionCoeffecient=1;
    mutationRate=0;
    carryingCapacity=1000;
    inbreedingCoeffecient=0;
    load=0;

    df <- data.frame(N=c(N),averageMutFreq=c(averageMutFreq),inbreedingCoeffecient=c(inbreedingCoeffecient))
    while (N<100) {
        intermediateAverageMutFreq = averageMutFreq - averageMutFreq^2 + inbreedingCoeffecient*averageMutFreq*(1-averageMutFreq);
        genedrift = rnorm(1,sd=sqrt(intermediateAverageMutFreq*(1-intermediateAverageMutFreq)/(2*mutations*N)));
        nAverageMutFreq = intermediateAverageMutFreq + genedrift;
        nInbreedingCoeffecient = (1-averageMutFreq)*(1/(2*N)+(1-1/(2*N))*inbreedingCoeffecient);
        # load = mutations*(averageMutFreq^2+inbreedingCoeffecient*averageMutFreq*(1-averageMutFreq));
        load = 1- (1-averageMutFreq^2)^mutations
        environmentStochasticity = 0;
        nN = N*exp(maxGrowthRate*(1-N/carryingCapacity)+environmentStochasticity)*(1-load) ;

        # t+1
        # append to dataframe
        N = nN;
        averageMutFreq = nAverageMutFreq;
        inbreedingCoeffecient = nInbreedingCoeffecient;
        df[nrow(df) + 1,] = c(N,averageMutFreq,inbreedingCoeffecient)
        if (is.nan(N)) {
            count=count+1
            break;
        }
        if (N<2) {
            count=count+1
            break;
        }
    }
}
print(1-count/replicates)
