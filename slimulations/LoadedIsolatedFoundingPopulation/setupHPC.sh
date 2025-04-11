outputFile=$1
if echo "$outputFile" | grep -q FitnessAverage; then
    echo "Seed,Fitness,GrowthRate,Time,Individuals,Females" >> "$outputFile"
else
    echo "Seed,Result,Time,Males,Individuals,GrowthRate,RecombinationRate,ChromosomeCount,MaxGenerations,CarryingCapacity,MutationFrequency,MutationCount,Sterile,Xlinked,FemaleOnlyEffect" >> "$outputFile"
fi


