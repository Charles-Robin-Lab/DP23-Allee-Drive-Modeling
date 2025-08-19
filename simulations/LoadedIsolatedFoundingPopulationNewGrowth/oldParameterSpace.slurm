#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# Source files
module load GCC/11.3.0
module load SLiM/4.3
source $bcpriSourcePath
# 10*11*11*7
# 15*15*30*14
# simulation count: 400*(4+(46-6)/8+1)*((102-2)/10+1)*((101-1)/10+1)*((0.31-0.01)/0.05+1)*2*2*2*1
# simulation count: 400*2*(4+2+(46-6)/8+1)*(10+(100-60)/10+1)*((300-10)/10+1)*((0.31-0.01)/0.05+1)*2*2*2
# simulation count: 400*(4+4+(40-10)/5+1)*(10+(100-60)/10+1)*((300-10)/10+1)*(10+4)*2*2*2
declare -A paramRanges=(
  ["carryingCapacity"]=1000
  ["chromosomeCount"]=300
  ["femaleOnlyEffect"]=1
  ["genomeSize"]="100e+6"
  ["maxGenerations"]=2000
  ["maximalAverageFemaleOffspring"]=$(seq 2.2 1.0 5.2; seq 6 8 46)
  ["migrantsSize"]=$(seq 2 10 102)
  ["mutsCount"]=$(seq 1 10 101)
  ["mutsFrequency"]=$(seq 0.01 0.05 0.31)
  ["postCompetitionMutationTiming"]=$(seq 0 1 1)
  ["recombinationRate"]=0
  ["sexed"]=1
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
)

paramOrdering=migrantsSize,maximalAverageFemaleOffspring,recombinationRate,chromosomeCount,maxGenerations,carryingCapacity,mutsFrequency,mutsCount,sterility,postCompetitionMutationTiming,xlinked,femaleOnlyEffect


runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulationNewGrowth.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

