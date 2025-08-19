#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath

declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=1
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]=3.0
  ["migrantsSize"]=25
  ["mutsCount"]=200
  ["mutsFrequency"]=0.04
  ["postCompetitionMutationTiming"]=0 
  ["recombinationRate"]="1e-3 5.62341325e-4 3.16227766e-4 1.77827941e-4 1e-4 5.62341325e-5 3.16227766e-5 1.77827941e-5 1e-5 5.62341325e-6 3.16227766e-6 1.77827941e-6 1e-6 5.62341325e-7 3.16227766e-7 1.77827941e-7 1e-7 5.62341325e-8 3.16227766e-8 1.77827941e-8 1e-8 1.77827941e-9 3.16227766e-9 5.62341325e-9 1e-9 5.62341325e-10 3.16227766e-10 1.77827941e-10 1e-10 5.62341325e-11 3.16227766e-11 1.77827941e-11 1e-11 5.62341325e-12 3.16227766e-12 1.77827941e-12 1e-12 5.62341325e-13 3.16227766e-13 1.77827941e-13 1e-13 5.62341325e-14 3.16227766e-14 1.77827941e-14 1e-14 5.62341325e-15 3.16227766e-15 1.77827941e-15 1e-15 0"
  ["sexed"]=1
  ["sterility"]="0 1"
  ["xlinked"]=0
  )

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

