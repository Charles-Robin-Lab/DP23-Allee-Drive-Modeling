#!/bin/bash 
# TODO:change this reference
source ../Bash-Command-Parameter-Ranges-Iterator/BashCommandParameterRangesIterator.sh

declare -A paramRanges=(
  ["migrantsSize"]=$(seq 2 1 30)
  ["mutsFrequency"]=$(seq 0.01 0.01 0.20)
  ["mutsCount"]=$(seq 3 1 30)
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
)

outputFile="out.csv"


# Check if the file already exists
if [ -f "$outputFile" ]; then
    echo "Error: File $outputFile already exists."
    exit 1
fi

# Write to the file
echo "Result,MutationFrequency,MutationCount,Individuals,Males,Sterile,Xlinked" >> "$outputFile"


runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
