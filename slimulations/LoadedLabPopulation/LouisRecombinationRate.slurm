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

# (4+4)*3*2*5*8*800
declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.5
  ["mutsCount"]=2
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=1
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.5
  ["mutsCount"]=2
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=2
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges


declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.33
  ["mutsCount"]=3
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=1
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.33
  ["mutsCount"]=3
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=3
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges


declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.25
  ["mutsCount"]=4
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=1
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.25
  ["mutsCount"]=4
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=4
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.2
  ["mutsCount"]=5
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=1
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

declare -A paramRanges=(
  ["maximalAverageFemaleOffspring"]=$(seq 2 1 5; seq 10 10 40)
  ["migrantsSize"]="2 6 10"
  ["mutsFrequency"]=0.2
  ["mutsCount"]=5
	["subpopsMutationsEvenlyDistributed"]=1
	["subpopsSexEvenlyDistributed"]=1
  ["carryingCapacity"]="100 1000"
  ["xlinked"]=0
  ["femaleOnlyEffect"]=0
  ["sterility"]=0
  ["chromosomeCount"]=5
  ["maxGenerations"]="5 10 20 100 2000"
  ["genomeSize"]=180000000
  ["recombinationRate"]="0"
)

runOverRanges "slim -d outputFilePath='$outputFile' %s LoadedLabPopulation.slim" "-d" paramRanges

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

