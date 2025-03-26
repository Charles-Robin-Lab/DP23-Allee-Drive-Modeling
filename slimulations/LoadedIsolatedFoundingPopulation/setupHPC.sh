outputFile=$1
if [["$outputFile" == *"FitnessAverage"*]]; then
    echo "Seed,Fitness,GrowthRate,Time,Individuals,Females"
else
    echo "Seed,Result,Time,Males,Individuals,GrowthRate,RecombinationRate,ChromosomeCount,MaxGenerations,CarryingCapacity,MutationFrequency,MutationCount,Sterile,Xlinked,FemaleOnlyEffect" >> "$outputFile"
fi


