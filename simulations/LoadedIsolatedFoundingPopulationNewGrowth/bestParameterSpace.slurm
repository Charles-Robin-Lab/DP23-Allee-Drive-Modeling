#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=20-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath
# simulation count: 400*((6-2.5)/0.5+1+3+3)*((120-30)/10+1+(25-5)/5+1)*((300-120)/20+1+(100-10)/10+1)*((0.15-0.01)/0.02+1)*2*2*2

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=300
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="$(seq 2.5 0.5 6) $(seq 10 5 20) $(seq 30 10 50)"
  ["migrantsSize"]="$(seq 5 5 25) $(seq 30 10 120)"
  ["mutsCount"]="$(seq 10 10 100) $(seq 120 20 300)"
  ["mutsFrequency"]="$(seq 0.01 0.01 0.05) $(seq 0.07 0.02 0.15)"
  ["postCompetitionMutationTiming"]=$(seq 0 1 1)
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
)

paramOrdering=,migrantsSize,maximalAverageFemaleOffspring,recombinationRate,chromosomeCount,maxGenerations,carryingCapacity,mutsFrequency,mutsCount,sterility,postCompetitionMutationTiming,xlinked,femaleOnlyEffect


runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges 

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

