# User Manual: Physiological Model for Algal 13C Fractionation

## Overview

This repository contains a mechanistic model of intracellular carbon fluxes and carbon isotope fractionation in algae. The model can be used to explore how carbon is exchanged among cellular compartments, how isotope fractionation develops during photosynthesis, and how environmental and physiological conditions influence model outputs.

The scripts are written for manual exploration and analysis. They can be used with either:
- synthetic or artificial input data for generic model experiments, or
- empirical culture data for comparison with observed fractionation patterns.

## Model structure

The model describes carbon transfer and isotope fractionation across four intracellular compartments:
1. Cytosol
2. Chloroplast stroma
3. Thylakoid
4. Pyrenoid

For haptophytes such as Gephyrocapsa (Emiliania huxleyi), calcification is treated as an additional process associated with the model framework rather than as a separate compartment.

## Repository files

The main files are:

- `R scripts/Carbon flux and isotope model.R`
  - Core model solver
  - Runs the carbon flux and isotope calculations for each simulation condition

- `R scripts/parameters.R`
  - Defines the model parameters and environmental inputs
  - Includes the values used for geometry, carbonate chemistry, pH, rate constants, and flux-related settings

- `R scripts/Generic model behavior.R`
  - Example script for exploring model behavior using artificial input data

- `R scripts/Model culture validation.R`
  - Example script for comparing model output with empirical culture data

## Required R packages

The scripts require the following R packages:

- `seacarb`
- `readxl`
- `ggplot2`
- `dplyr`

In addition, the ODE solver file `grind.R` (De Boer, 2024) (https://bioinformatics.bio.uu.nl/rdb/grind.html) is required. It should be placed in the working directory used to run the scripts.

## Working directory setup

A convenient way to run the model is to place the `data` folder inside the `R scripts` folder and set the working directory to that folder. This makes it easier to reference data files using relative paths such as:

```r
source("Generic model behavior.R")
```

or

```r
source("Model culture validation.R")
```

This setup is especially helpful when working locally or when sharing the repository with others.

## Running the model

### 1. Explore generic model behavior

Use the generic behavior script when you want to investigate how the model responds to changes in environmental or physiological conditions.

```r
source("Generic model behavior.R")
```

This script creates an example data frame called `yourdata` and runs the model across a set of conditions. The `yourdata` object in this script is only a simple dummy example intended to show the workflow. You can replace it with any custom input table that matches the model requirements, and you can vary any parameter defined in `parameters.R` to explore how the model behaves under different assumptions.

In other words, the example data frame is not meant to be a fixed or restrictive setup. It is simply a starting point for exploring model behavior. You can change:
- environmental inputs such as CO2, light, temperature, salinity, and pH
- physiological inputs such as cell volume, growth rate, and carbon fixation rate
- any parameter in `parameters.R`, including geometry, permeability, pH values, reaction rates, and carbonate chemistry settings

### 2. Compare the model with culture data

Use the culture validation script when you want to compare the model with measured fractionation values.

```r
source("Model culture validation.R")
```

The script reads culture data from the `data` folder and runs the model for each condition. It then compares model outputs with observed values and produces plots for inspection.

## Parameter exploration

A key feature of this model is that the parameters defined in `R scripts/parameters.R` can be changed freely to investigate model behavior.

This includes:
- Rubisco fractionation values
- membrane permeability
- cell geometry
- compartment pH values
- reaction rate constants
- carbonate chemistry settings
- carbon fixation rate inputs
- uptake-related and leakage-related parameters

In other words, the parameter file is not limited to a small set of values that are treated as “tunable” for fitting. The non-tunable values can also be edited manually to test how assumptions about physiology and chemistry influence the model outputs.

This makes the model suitable for:
- sensitivity analysis
- conceptual experiments
- manual testing of alternative parameter assumptions
- exploration of how different biological processes affect isotope fractionation

## Input data format

The scripts use a data frame containing one row per simulation condition. The data frame can be built manually or loaded from Excel files in the `data` folder.

Typical columns include:
- `daylight`
- `co2`
- `pH`
- `Vcell`
- `eps`
- `Temp`
- `sal`
- `PFD`
- `species`

Additional columns may be used depending on the specific script and the information available for the run.

## Output variables

The model returns outputs such as:
- photosynthetic fractionation (`eps`)
- calcification-related fractionation (`EpcalcCO2`)
- isotope values for carbon species in different compartments
- relative HCO3- uptake
- CO2 leakage
- external CO2 concentration and related variables

These outputs can be used to interpret how changing environmental or physiological assumptions affects the model.

## Troubleshooting

Common issues and suggested fixes:

- Package installation errors
  - Install the required packages before running the scripts:

```r
install.packages(c("seacarb", "readxl", "ggplot2", "dplyr"))
```

- `grind.R` not found
  - Make sure `grind.R` is present in the working directory or in a directory that is on your R search path.
  - If you are using RStudio, confirm that the correct folder is selected as the working directory.

- File not found errors for input files
  - Check that the `data` folder exists and contains the expected files.
  - Verify that your working directory is set to the folder containing the scripts or to the folder that contains `data`.
  - You can inspect the current working directory with:

```r
getwd()
```

- Errors when sourcing scripts
  - Ensure that you run the script from the correct directory, or use `setwd()` to point R to the appropriate folder before sourcing:

```r
setwd("path/to/your/R scripts/folder")
source("Generic model behavior.R")
```

- Unexpected or empty model output
  - Check that the input data frame contains the expected columns and that the values are valid.
  - Review the parameter values in `R scripts/parameters.R` for obvious inconsistencies such as missing values or unrealistic ranges.
  - If the script stops with an error, read the message carefully because it often indicates the first missing object, package, or file.

## Notes

- The scripts are intended for manual exploration and analysis.
- The model can be run with either synthetic input data or empirical culture data.
- The parameter file is the main place to modify the model setup for new experiments.
- For convenience, the `data` folder may be placed inside the `R scripts` folder so that relative file paths are simpler to manage.
