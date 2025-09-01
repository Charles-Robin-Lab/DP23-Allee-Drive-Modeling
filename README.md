# Genetic Allee effects for controlling invasive founder populations
 This repository models Genetic Allee Effects in founder populations arising from a core population with a given load. See the full manuscript here: [10.22541/au.172114604.45747531/v1](https://doi.org/10.22541/au.172114604.45747531/v1)

## Run Jobs/Generate Data
These jobs are set up to run on the University of Melbourne's High Performance Computing system "Spartan". You will need to modify `simulations/slimSBatch` if you are using a different SLURM environment.

To run a job:
```sh
cd simulations
./slimSBatch.sh LoadedIsolatedFoundingPopulationNewGrowth/[JobName].slurm
```

## Run data analysis/Generate base figures
Put the appropriate data from [here](https://drive.google.com/drive/folders/1GEInfEuVN_sUtVFbAjZL5v0fcvSI2w4c?usp=drive_link) into the `data/` directory. R code in the `data_analysis/` directory will produce the base figures for the paper. This code should be run with the working directory set to the top of the repo. Code for specific figures can be found by grepping for `figures/` i.e.
```sh
grep -r "figures/*" data_analysis/
```

## Assemble figures
SVG output from the data_analysis code and other printed output is used to assemble the figures. Affinity designer and GIMP files used during assembly can be found in the `/figures` directory. In some cases the base SVG will need to be edited to fix imcompadibility between affinity designer and svglite.
