#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)
# Edited heavily by Louis Nowell Nicolle
# This file will generate the data used for figures: 4

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath
# simulation count: 400*2*4*((0.10-0.0)/0.0025+1)*((400-80)/10+1)
declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=400
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]="3 5"
  ["migrantsSize"]="12 25 50 100"
  ["mutsCount"]=$(seq 80 10 400)
  ["mutsFrequency"]=$(seq 0.00 0.0025 0.10)
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=0
  ["xlinked"]=0
)

paramOrdering=,migrantsSize,maximalAverageFemaleOffspring,recombinationRate,chromosomeCount,maxGenerations,carryingCapacity,mutsFrequency,mutsCount,sterility,postCompetitionMutationTiming,xlinked,femaleOnlyEffect
runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges 
##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

