#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-800

# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.0.1
source $bcpriSourcePath

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=3.0
  ["migrantsSize"]=10
  ["mutsFrequency"]=0.05
  ["mutsCount"]=40
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]=100
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=2
  ["genomeSize"]=180000000
  ["recombinationRate"]="2.32e-8 4.22e-8"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

