#!/bin/bash 
# TODO:change this reference
source ../../Coding/Bash-Command-Parameter-Ranges-Iterator/BashCommandParameterRangesIterator.sh

declare -A paramRanges=(
  ["migrantsSize"]=$(seq 2 5 30),
  ["mutsFrequency"]=$(seq 0.01 0.05 0.20),
  ["mutsCount"]=$(seq 3 4 30),
  ["sterility"]=$(seq 0 1 1),
  ["xlinked"]=$(seq 0 1 1)
)

runOverRanges "slim %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges