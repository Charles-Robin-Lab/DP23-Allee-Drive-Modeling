#!/usr/local/easybuild-2019/easybuild/software/compiler/gcc/11.2.0/bash/4.3/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="LoadedIsolatedFoundingPopulation"

# The project ID which this job should run under:
#SBATCH --account="punim2001"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Use this email address:
#SBATCH --mail-user=louis.nowellnicolle@student.unimelb.edu.au

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# ends successfully
#SBATCH --mail-type=END

# log output files
#SBATCH -o ../logs/slurm.%N.%j.out # STDOUT 
#SBATCH -e ../logs/slurm.%N.%j.err # STDERR

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-400

# # # File and version checks 

# Run the job from this directory:
cd /home/lnowellnicol/punim2001/DP23-Allee-Drive-Modeling

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Check if the bash version is lower than 4.3.0 since we use associative arrays
if [[ "${BASH_VERSION}" < "4.3.0" ]]; then
    echo "Error: Bash version 4.3.0 or higher is required. Please edit the first line of the script to point to it"
    exit 1
fi

# Write csv header
outputFile="out_$SLURM_ARRAY_JOB_ID.csv"
if [ "$SLURM_ARRAY_TASK_ID" == "1" ]; then
   echo "Seed,Result,Time,MutationFrequency,MutationCount,GrowthRate,RecombinationRate,Individuals,Males,Sterile,Xlinked" >> "$outputFile"
fi
# Wait for csv to exist
while [ ! -f "$outputFile" ]; do
  :
done

sleep 1

# check sources
bcpriSourcePath="../Bash-Command-Parameter-Ranges-Iterator/BashCommandParameterRangesIterator.sh"
if [ ! -f "$bcpriSourcePath" ]; then
    echo "Error: Bash-Command-Parameter-Ranges-Iterator not found please download from https://github.com/4321louis/Bash-Command-Parameter-Ranges-Iterator and edit this script to point to it"
    exit 1
fi

# # # run job code

# Source files
module load slim/4.0.1
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

