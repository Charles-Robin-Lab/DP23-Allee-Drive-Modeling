#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)
# Edited heavily by Louis Nowell Nicolle
# This file will generate the data used for figures: 3, S3

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath

# 85% fitness in core
declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=200
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=50
  ["maximalAverageFemaleOffspring"]="3"
  ["migrantsSize"]=$(seq 5 1 75)
  ["mutsCount"]=200
  ["mutsFrequency"]=0.04
  ["postCompetitionMutationTiming"]="0 1"
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]="0 1"
  ["xlinked"]=0
)
runOverRanges "slim -d fitnessLogFilePath='${outputFile}' -d outputFilePath='${outputFile}_params.csv' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges
##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

