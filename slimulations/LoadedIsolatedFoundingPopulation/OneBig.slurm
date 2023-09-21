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
module load SLiM/4.0.1
source $bcpriSourcePath

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=3.0
  ["migrantsSize"]=25
  ["mutsFrequency"]=0.05
  ["mutsCount"]=100
  ["xlinked"]=0
  ["sterility"]=1
  ["chromosomeCount"]=1
  ["genomeSize"]=100000000
  ["recombinationRate"]="1e-6 1.77827941e-6 3.16227766e-6 5.62341325e-6 1e-7 1.77827941e-7 3.16227766e-7 5.62341325e-7 1e-8 1.77827941e-8 3.16227766e-8 5.62341325e-8 1e-9 1.77827941e-9 3.16227766e-9 5.62341325e-9 1e-10 1.77827941e-10 3.16227766e-10 5.62341325e-10 1e-11 1.77827941e-11 3.16227766e-11 5.62341325e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s slimulations/LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

