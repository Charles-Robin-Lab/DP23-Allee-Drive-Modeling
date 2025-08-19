#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

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
  ["chromosomeCount"]=250
  ["femaleOnlyEffect"]=0
  ["genomeSize"]=1000
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]=$(seq 1.0 0.25 3.0)
  ["migrantsSize"]=$(seq 10 10 100)
  ["mutsCount"]=$(seq 10 40 250)
  ["mutsFrequency"]=$(seq 0.01 0.01 0.15)
  ["postCompetitionMutationTiming"]=0
  ["recombinationRate"]=0
  ["sexed"]=0
  ["sterility"]=0
  ["xlinked"]=0
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

