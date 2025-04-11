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
#SBATCH -o logs/slurm.%N.%j.out # STDOUT 
#SBATCH -e logs/slurm.%N.%j.err # STDERR

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
# Parralelise the computation
#SBATCH --array=1-100

# # # File and version checks 

# Run the job from this directory:
cd /home/lnowellnicol/punim2001/DP23-Allee-Drive-Modeling




echo "$SLURM_JOB_ID $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $SLURM_JOB_START_TIME $((SLURM_ARRAY_JOB_ID-SLURM_ARRAY_TASK_ID))">> test

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

