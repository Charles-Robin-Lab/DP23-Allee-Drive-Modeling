#! /bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu May 18 2023 02:25:29 GMT+1000 (Australian Eastern Standard Time)
# Edited by Louis Nowell Nicolle
# This file will generate the data used for figures: S2

# # # Settings for Sbatch

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-8:00:00
# Parralelise the computation
#SBATCH --array=1-400


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
#SBATCH -o /dev/null
#SBATCH -e /home/lnowellnicol/punim2001/DP23-Allee-Drive-Modeling/logs/slurm.numericalSolutionPoint1.err

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load foss/2022a R/4.2.2

# Run the job from the directory where it was launched (default)

# The job command(s):
Rscript ./slimulations/Numerical/bestModelDiffs.r bestmodeldiff_02072025_2

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s