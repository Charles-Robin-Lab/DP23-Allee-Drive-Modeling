#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-1000

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.0.1
source $bcpriSourcePath

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=3
  ["migrantsSize"]=25
  ["mutsFrequency"]=0.075
  ["mutsCount"]=80
  ["sterility"]=0
)
runOverRanges "slim -d fitnessLogFilePath='$outputFile' outputFilePath='/dev/null' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=3.5
  ["migrantsSize"]=25
  ["mutsFrequency"]=0.075
  ["mutsCount"]=80
  ["sterility"]=0
)

runOverRanges "slim -d fitnessLogFilePath='$outputFile' outputFilePath='/dev/null'  %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=4
  ["migrantsSize"]=25
  ["mutsFrequency"]=0.075
  ["mutsCount"]=80
  ["sterility"]=0
)

runOverRanges "slim -d fitnessLogFilePath='$outputFile' outputFilePath='/dev/null'  %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

