#!/bin/bash

# What file is being run?
slurmscript=$(realpath $1)
if [ ! -n $slurmscript ]; then
echo "please provide full path to a slurm job script as argument" >&2
exit 1
fi
situationName="$(echo $slurmscript | sed -r "s/.+\/(.+)\..+/\1/" )"
modelName="$(dirname $slurmscript | sed 's:.*/::' )"
jobName="${situationName}_$(modelName | sed 's/[a-z]//g')_$(date +%s)"


# Partition for the job:
sbatchOptions="--partition=cascade"

# Multithreaded (SMP) job: must run on one node 
sbatchOptions+=" --nodes=1"

# The project ID which this job should run under:
sbatchOptions+=" --account=\"punim2001\""

# Maximum number of tasks/CPU cores used by the job:
sbatchOptions+=" --ntasks=1"
sbatchOptions+=" --cpus-per-task=1"

# Use this email address:
sbatchOptions+=" --mail-user=louis.nowellnicolle@student.unimelb.edu.au"

# Send yourself an email when the job:
# aborts abnormally (fails)
sbatchOptions+=" --mail-type=FAIL"
# ends successfully
sbatchOptions+=" --mail-type=END"

# log output files
sbatchOptions+=" -o ../logs/slurm.$jobName.out" # STDOUT 
sbatchOptions+=" -e ../logs/slurm.$jobName.err" # STDERR


# # # File and version checks 
# Run the job from this directory:
cd /home/lnowellnicol/punim2001/DP23-Allee-Drive-Modeling

# Check if the bash version is lower than 4.3.0 since we use associative arrays
if [[ "${BASH_VERSION}" < "4.3.0" ]]; then
    echo "Error: Bash version 4.3.0 or higher is required. Please edit the first line of the script to point to it">&2
    exit 1
fi

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


# Write csv header
outputFile="./data/out_${jobName}.csv"
echo "Seed,Result,Time,MutationFrequency,MutationCount,GrowthRate,RecombinationRate,Individuals,Males,Sterile,Xlinked" >> "$outputFile"

# Runscript
workdir=$(dirname $slurmscript)
cmd="sbatch ${sbatchOptions} --chdir=$workdir $slurmscript"
echo $cmd
$cmd

