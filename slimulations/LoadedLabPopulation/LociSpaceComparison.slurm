#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-0:0:00
# Parralelise the computation
#SBATCH --array=1-200

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.0.1
source $bcpriSourcePath

#is the distribution of the locis the cause/ turn off the distribution toggle


# 800*(4+2)*7*8*(4+7)*2
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 20)
  ["migrantsSize"]=$(seq 10 15 100)
  ["mutsFrequency"]=$(seq 0.05 0.05 0.4)
  ["mutsCount"]=$(seq 2 1 5; seq 10 15 100)
	["subpopsMutationsEvenlyDistributed"]=0
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=1
  ["maxGenerations"]="2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 20)
  ["migrantsSize"]=$(seq 10 15 100)
  ["mutsFrequency"]=$(seq 0.05 0.05 0.4)
  ["mutsCount"]=$(seq 2 1 5; seq 10 15 100)
	["subpopsMutationsEvenlyDistributed"]=0
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=100
  ["maxGenerations"]="2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

