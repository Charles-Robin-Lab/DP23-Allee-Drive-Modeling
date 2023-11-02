#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-0:0:00
# Parralelise the computation
#SBATCH --array=1-800
replicates=50
# # # run job code

# Source files
module load GCC/11.3.0
module load SLiM/4.0.1
source $bcpriSourcePath
#50*(4+4)*5*6*7*4*16 
RR=('0' '1e-5' '2.32e-8' '4.22e-8' '1e-8' 1e-12)

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]=$(seq 2 2 10)
  ["mutsFrequency"]=$(seq 0.05 0.05 0.30)
  ["mutsCount"]=$(seq 3 3 21)
  ["sterility"]=0
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["genomeSize"]=180e+6
  ["carryingCapacity"]=100
  ["chromosomeCount"]=3
  ["maxGenerations"]="5 10 20 100 2000"
  ["recombinationRate"]=${RR[$SLURM_ARRAY_TASK_ID % 16]}
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges
##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

