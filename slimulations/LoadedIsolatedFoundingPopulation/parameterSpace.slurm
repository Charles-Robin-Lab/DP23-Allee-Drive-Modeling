#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=3-0:0:00
# Parralelise the computation
#SBATCH --array=1-1200

# Source files
module load SLiM/4.0.1
source $bcpriSourcePath
# full
# declare -A paramRanges=(
#   ["maximalGrowthRate"]=$(seq 0.1 0.1 1.5) $(seq 2 0.5 10)
#   ["migrantsSize"]=$(seq 2 1 40)
#   ["mutsFrequency"]=$(seq 0.01 0.01 0.20)
#   ["mutsCount"]=$(seq 3 1 100)
#   ["sterility"]=$(seq 0 1 1)
#   ["xlinked"]=$(seq 0 1 1)
# )
# intial test
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2.2 1.0 5.2; seq 6 8 46)
  ["migrantsSize"]=$(seq 2 10 102)
  ["mutsFrequency"]=$(seq 0.01 0.05 0.31)
  ["mutsCount"]=$(seq 1 10 101)
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

