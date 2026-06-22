# Physiological Model for Algal 13C Fractionation

A mechanistic model of intracellular carbon fluxes and carbon isotope fractionation in algae, simulating the complex pathways of carbon through different cellular compartments. The model can be used generically to investigate model behavior or fitted against empirical culture data.

## Overview

This repository contains a comprehensive model of carbon isotope fractionation in algal photosynthesis, with applications to:
- **Gephyrocapsa (Emiliania) huxleyi** - A coccolithophore species included with culture validation data
- **Dinoflagellates** - Including *Alexandrium tamarense*, *Scrippsiella trochoidea*, *Gonyaulax spinifera*, and *Protoceratium reticulatum*

The model simulates:
1. **Carbon flux dynamics** through cellular compartments (cytosol, chloroplast, thylakoid, pyrenoid)
2. **Carbon isotope fractionation** (δ¹³C) during photosynthetic carbon fixation and transport
3. **Response to environmental conditions** including CO₂ concentration, light intensity, temperature, salinity, and pH

## Repository Structure

### Scripts

- **`Carbon flux and isotope model.R`** - Core model solver
  - Implements ODE system for carbon concentrations and isotope ratios
  - Calculates compartment-specific δ¹³C values
  - Determines photosynthetic fractionation (ε_p) and calcification fractionation (ε_PIC)
  - Returns model outputs for individual culture conditions

- **`parameters.R`** - Parameter configuration file
  - Rubisco fractionation values (species-specific)
  - Membrane permeability coefficients
  - Cell geometry calculations (cytosol, chloroplast, pyrenoid, thylakoid dimensions)
  - Reaction rate constants for CO₂/HCO₃⁻ inter-conversion
  - Carbonate chemistry calculations
  - Tunable parameters (calcification rate, uptake rates, Michaelis-Menten parameterization)

- **`Generic model behavior.R`** - Model investigation script
  - Explores model behavior across environmental gradients (example demonstrates CO₂ and light intensity, but any parameter can be varied)
  - Creates artificial data to study model sensitivity and response patterns
  - Generates visualization of modeled outputs (e.g., photosynthetic fractionation)
  - Useful for understanding model behavior without empirical data constraints

- **`Model culture validation.R`** - Model-data comparison script
  - Loads empirical culture data (either *E. huxleyi* or dinoflagellates)
  - Compares model outputs to measured fractionation
  - Calculates goodness-of-fit metrics (RMSE, R²)
  - **Post hoc validation plots on independent culture data:**
    - Photosynthetic fractionation (ε_p) vs CO₂
    - HCO₃⁻ uptake fraction vs CO₂
    - CO₂ leakage vs CO₂
    - Calcification fractionation (ε_PIC) vs light intensity

### Data Files

**E. huxleyi culture data:**
- `Ehux.xlsx` - Culture measurements of δ¹³C (bulk and PIC), growth rate, cell geometry, and environmental conditions
- `haptfluxdata.xlsx` - Carbon flux measurements (HCO₃⁻ uptake fraction and CO₂ leakage)

**Dinoflagellate culture data:**
- `Dino.xlsx` - Culture measurements for *A. tamarense*, *S. trochoidea*, *G. spinifera*, and *P. reticulatum*
- `dinofluxdata.xlsx` - Carbon flux measurements for dinoflagellate species

## Required Packages

```R
seacarb          # Carbonate chemistry calculations
readxl           # Reading Excel data files
ggplot2          # Data visualization
dplyr            # Data manipulation
grind.R*         # ODE solver (external file)
```

*Download `grind.R` from: https://bioinformatics.bio.uu.nl/rdb/grind.html

## Quick Start

### 1. Investigate Generic Model Behavior

```R
source("Generic model behavior.R")
```

This script:
- Creates artificial parameter gradients (by default CO₂ and light intensity, but customizable)
- Runs the model across these conditions
- Visualizes predicted model outputs

### 2. Validate Against Culture Data

```R
# For E. huxleyi:
source("Model culture validation.R")

# For dinoflagellates (modify line 48 in Model culture validation.R):
# culture.data <- as.data.frame(read_excel("data/Dino.xlsx"))
# fluxdata <- read_excel("data/dinofluxdata.xlsx")
```

This script performs post hoc validation by:
- Loading independent culture data
- Running the model with culture-specific parameters
- Comparing modeled vs measured fractionation
- Evaluating model fit (RMSE, R²)
- Plotting metabolic fluxes (HCO₃⁻ uptake, CO₂ leakage) against independent measurements

## Data Format Requirements

Culture data files should contain these columns:
- `daylight` - Daylength (h)
- `co2` - CO₂ concentration (µmol/L)
- `pH` - pH
- `Vcell` - Cell volume (cm³)
- `eps` - Photosynthetic fractionation ε_p (‰)
- `Temp` - Temperature (°C)
- `sal` - Salinity (PSU)
- `PFD` - Photon flux density (µmol photons/m²/s)
- `POC` - Cellular carbon content (mol) *if not providing r1*
- `ui` - Instantaneous growth rate (/day) *if not providing r1*
- `r1` - Carbon fixation rate (mol/cm³/s) *if not providing POC and ui*
- `species` - Species name (must match parameter settings)
- `EpcalcCO2` - Calcification fractionation (‰)
- `ref` - (optional) Reference for data source

## Model Structure

### Compartments

The model simulates carbon dynamics across five intracellular compartments:
1. **Cytosol** - Outer cellular compartment, site of initial CO₂/HCO₃⁻ exchange
2. **Chloroplast stroma** - Inner chloroplast compartment, connected to thylakoids
3. **Thylakoid** - Site of photosynthetic light reactions/ATP production
4. **Pyrenoid** - RuBisCO-containing compartment, site of carbon fixation
5. **Calcification** - CaCO₃ precipitation (for *E. huxleyi*)

### Metabolic Fluxes

Key processes modeled:
- CO₂/HCO₃⁻ inter-conversion at each compartment boundary
- CO₂ diffusion across membranes
- HCO₃⁻ active uptake (saturable Michaelis-Menten kinetics)
- Carbon fixation by RuBisCO
- CaCO₃ formation (coccolithophore-specific)

### Fractionation Mechanisms

The model accounts for isotopic fractionation during:
1. **Enzymatic fixation** - RuBisCO fractionation (~11‰ for *E. huxleyi*, ~24‰ for dinoflagellates)
2. **CO₂/HCO₃⁻ equilibrium** - Temperature-dependent isotopic equilibrium
3. **Transport** - Diffusive and active transport across membranes
4. **Calcification** - CaCO₃ precipitation fractionation

## Model Parameters

### Fixed Parameters

| Parameter | E. huxleyi | Dinoflagellates | Description |
|-----------|-----------|-----------------|-------------|
| RuBisCO fractionation | 11.1‰ | 24‰ | Isotopic fractionation during CO₂ fixation |
| Membrane permeability | 0.01 cm/s | 0.027 cm/s | CO₂ diffusion across biological membranes |
| Vcell | | | Constrained from culture measurements|
| pH | | | Constrained from culture measurements|
| Temp and sal | | | Constrained from culture measurements|
| r1/POC/ui | | | growth rate --> constrained from culture measurements|
| PFD | | | Constrained from culture measurements|
| Daylight | | | Constrained from culture measurements|






### Tunable Parameters

These are optimized to fit culture data:
- `calc` - Calcification rate (*E. huxleyi* only)
- `uptake.rate` - Active HCO₃⁻ uptake rate (species-specific)
- `Kn` - Michaelis-Menten constant for HCO₃⁻ uptake (*E. huxleyi* only)
- `pP` - Pyrenoid permeability factor

## Output Variables

The model returns a data frame with these columns for each simulation:

| Variable | Unit | Description |
|----------|------|-------------|
| `eps` | ‰ | Photosynthetic fractionation (ε_p) |
| `EpcalcCO2` | ‰ | Calcification fractionation (ε_PIC) |
| `dCc`, `dHc` | ‰ | δ¹³C of CO₂ and HCO₃⁻ in cytosol |
| `dCch`, `dHch` | ‰ | δ¹³C in chloroplast stroma |
| `dCt`, `dHt` | ‰ | δ¹³C in thylakoid |
| `dCp`, `dHp` | ‰ | δ¹³C in pyrenoid |
| `Cc`, `Cch`, `Ct`, `Cp` | mol/cm³ | Carbon concentrations in compartments |
| `Uptake` | fraction | Relative HCO₃⁻ uptake (0-1) |
| `Leakage` | fraction | CO₂ leakage rate (0-1) |
| `Vcell` | cm³ | Cell volume |
| `Ce` | mol/cm³ | External CO₂ concentration |
| `pCO2` | µatm | Partial pressure of CO₂ |
| `Temp` | °C | Temperature |
| `PFD` | µmol/m²/s | Photon flux density |

## License

MIT License - See LICENSE file for details

## Citation

If you use this model in your research, please cite the repository and any associated publications.

## Contact

For questions or issues, please open an issue on this repository.
