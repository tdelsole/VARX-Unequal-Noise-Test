# VARX-Unequal-Noise-Test
This repository contains the R codes and data used in the paper:   **"Testing Equality of Autoregressive Parameters Without Assuming Equality of Noise Variances"**   by *Timothy DelSole* and *Michael K. Tippett*.  
The paper proposes a **significance test** for determining if the parameters of two VARX (Vector AutoRegressive model with eXogenous inputs) models differ, **without assuming equality of noise covariances**.

---

## Data

- `ipcc.global.ts_Amon_CMIP6_historical.consolidate.RData`:  
  Monthly 2m-temperature from **CMIP6 historical runs**, projected onto the 58 Atlas regions defined in [Iturbide et al., 2020].

- `ipcc.global.2t_Amon_CMIP6_ERA5.consolidate2.RData`:  
  Monthly 2m-temperature from **ERA5 reanalysis** ([Hersbach et al., 2020]), also projected onto the 58 Atlas regions.

- `IPCC-WGI-reference-regions-v4_R.rda`:  
  Edges of polygons for the 58 Atlas regions defined in [Iturbide et al., 2020].
---

## Workflow

To reproduce the analysis:  

1. Run `NoiseNotFirst.Github.R`.  
2. Edit **line 15** to set one of the two cases:  
   - `first.say = 'ar'`  
   - `first.say = 'forcing'`  
3. The script:
   - Reads the data.
   - Computes the 2m-T averaged over 5 regions.
   - Compares VARX models **without assuming equality of noise covariances**.  
4. Results are saved as figures in the `./figures/` directory.

---

## Code Overview

| File | Description |
|------|-------------|
| `NoiseNotFirst.Github.R` | Main script; reads data and runs tests for every CMIP6 historical simulation. |
| `mic.master.R` | Computes the Mutual Information Criterion (MIC) to select model order. |
| `anncyc.remove.monthly.R` | Removes the annual cycle from a time series array. |
| `pdf.eps.R` | Produces PDF output files. |
| `plot_latlon_v4.R` | Creates spatial plots. |
| `AtlasMaskProduce.R` | Defines the 58 Atlas regions for CMIP6. |
| `diff.regression.nested.mult.NoiseNotFirst.R` | Core function for testing equality of VARX parameters. |
| `deviance.AnyAR.FixNoise.NoiseNotFirst.R` | High-level function for data formatting and hypothesis testing. |
| `gev.R` | Solves the generalized eigenvalue problem. |
| `LjungBox.R` | Performs the multivariate Ljung-Box test. |
| `n.to.monthly.R` | Maps one set of date indices to another. |
| `simulate.AnyAR.FixNoise.1forcing.R` | Simulates a VARX model from estimated parameters. |
| `sync.periods.R` | Synchronizes datasets with overlapping time domains. |
| `timeseries2VARX.cyclo.R` | Formats time series data for regression analysis. |

---

## Quick Start

```r
# Run the main script from the R console
source("NoiseNotFirst.Github.R")

---

## Citation

If you use these codes or data, please cite:
DelSole, T. and Tippett, M.K. (2025). Testing Equality of Autoregressive Parameters Without Assuming Equality of Noise Variances.


