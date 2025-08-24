#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)
# Edited heavily by Louis Nowell Nicolle
# This file will generate the data used for figures: S10

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=3-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath



declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]="1 100"
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="3 3.5 4"
  ["migrantsSize"]=$(seq 10 5 100)
  ["mutsCount"]=100
  ["mutsFrequency"]=0.083
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=0
  ["xlinked"]=0
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]="1 100"
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="3 3.5 4"
  ["migrantsSize"]=$(seq 10 5 100)
  ["mutsCount"]=70
  ["mutsFrequency"]=0.099
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=0
  ["xlinked"]=0
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]="1 100"
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="3 3.5 4"
  ["migrantsSize"]=$(seq 10 5 100)
  ["mutsCount"]=40
  ["mutsFrequency"]=0.131
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=0
  ["xlinked"]=0
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]="1 100"
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="3 3.5 4"
  ["migrantsSize"]=$(seq 10 5 100)
  ["mutsCount"]=10
  ["mutsFrequency"]=0.259
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=0
  ["xlinked"]=0
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

