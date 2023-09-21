#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)

# # # Settings for Sbatch

# Partition for the job:
#SBATCH --partition=cascade

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
#SBATCH --array=1-800

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
outputFile="./data/out_FourNormal400_$SLURM_ARRAY_JOB_ID.csv"
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
module load SLiM/4.0.1
source $bcpriSourcePath

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=3.0
  ["migrantsSize"]=25
  ["mutsFrequency"]=0.015
  ["mutsCount"]=400
  ["xlinked"]=0
  ["sterility"]=1
  ["chromosomeCount"]=4
  ["genomeSize"]=100000000
  ["recombinationRate"]="1e-6 1.77827941e-6 3.16227766e-6 5.62341325e-6 1e-7 1.77827941e-7 3.16227766e-7 5.62341325e-7 1e-8 1.77827941e-8 3.16227766e-8 5.62341325e-8 1e-9 1.77827941e-9 3.16227766e-9 5.62341325e-9 1e-10 1.77827941e-10 3.16227766e-10 5.62341325e-10 1e-11 1.77827941e-11 3.16227766e-11 5.62341325e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s slimulations/LoadedIsolatedFoundingPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

