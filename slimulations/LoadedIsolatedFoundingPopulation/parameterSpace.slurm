#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=3-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# Source files
module load GCC/11.3.0
module load SLiM/4.0.1
source $bcpriSourcePath
# simiulation count: 400*(4+(46-6)/8+1)*((102-2)/10+1)*((0.31-0.01)/0.05+1)*((101-1)/10+1)*2*2*1
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2.2 1.0 5.2; seq 6 8 46)
  ["migrantsSize"]=$(seq 2 10 102)
  ["mutsFrequency"]=$(seq 0.01 0.05 0.31)
  ["mutsCount"]=$(seq 1 10 101)
  ["sterility"]=$(seq 0 1 1)
  ["xlinked"]=$(seq 0 1 1)
  ["recombinationRate"]="1e-8"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

