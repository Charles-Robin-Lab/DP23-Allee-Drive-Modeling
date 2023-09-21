#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# # # run job code

# Source files
module load GCC/11.3.0
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
# declare -A paramRanges=(
#  ["maximalGrowthRate"]=$(seq 0.1 0.5 1.6; seq 2 4 22)
#  ["migrantsSize"]=$(seq 2 10 102)
#  ["mutsFrequency"]=$(seq 0.01 0.05 0.31)
#  ["mutsCount"]=$(seq 1 10 101)
#  ["sterility"]=$(seq 0 1 1)
#  ["xlinked"]=$(seq 0 1 1)
#)


# nice graphs
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2.2 0.04 7)
  ["migrantsSize"]=20
  ["mutsFrequency"]=0.1
  ["mutsCount"]=60
  ["sterility"]=$(seq 0 1 1)
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=4.0
  ["migrantsSize"]=$(seq 2 1 120)
  ["mutsFrequency"]=0.1
  ["mutsCount"]=60
  ["sterility"]=$(seq 0 1 1)
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=4.0
  ["migrantsSize"]=20
  ["mutsFrequency"]=$(seq 0.04 0.002 0.15)
  ["mutsCount"]=60
  ["sterility"]=$(seq 0 1 1)
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=4.0
  ["migrantsSize"]=20
  ["mutsFrequency"]=0.1
  ["mutsCount"]=$(seq 1 1 120)
  ["sterility"]=$(seq 0 1 1)
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

