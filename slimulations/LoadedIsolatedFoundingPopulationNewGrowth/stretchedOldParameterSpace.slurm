#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath
# simulation count: 400*(5+(46-6)/10+1)*((122-2)/10+1)*((290-10)/20+1)*((0.13-0.01)/0.03+1)*2*2*2

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=300
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="2.2 2.95 3.7 4.45 5.2 $(seq 6 10 46)"
  ["migrantsSize"]=$(seq 2 10 122)
  ["mutsCount"]=$(seq 10 20 290)
  ["mutsFrequency"]=$(seq 0.01 0.03 0.13)
  ["postCompetitionMutationTiming"]=$(seq 0 1 1)
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
)

paramOrdering=,migrantsSize,maximalAverageFemaleOffspring,recombinationRate,chromosomeCount,maxGenerations,carryingCapacity,mutsFrequency,mutsCount,sterility,postCompetitionMutationTiming,xlinked,femaleOnlyEffect


runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges $paramOrdering /home/lnowellnicol/punim2001/DP23-Allee-Drive-Modeling/data/out_stretchedOldParameterSpace_LIFPNG_1750521215.csv 400

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

