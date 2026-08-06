# Correcting Regression to the Mean

This repository contains the code, data, and generated figures for the
manuscript "Failure to distinguish between measurement error and change
heterogeneity leads to improper correction of regression-to-the-mean artifacts,"
accepted for publication in *The American Naturalist*.

The analyses evaluate when apparent relationships between trait change and
initial trait values arise from regression-to-the-mean artifacts, and compare
several approaches for estimating or correcting those relationships.

## Repository Manifest

### Main

| File | Description |
|---|---|
| `README.md` | Repository documentation and project overview |
| `RegToMean.Rproj` | RStudio project file |

### Code

| File | Description |
|---|---|
| `code/RTM_analysis.Rmd` | Canonical analysis notebook for simulations, empirical analyses, and manuscript figure generation; defines the helper functions used directly in the main analysis |
| `code/acer_variance.R` | Calculates measurement-error variance for the *A. cervicornis* case study; this value is carried into `RTM_analysis.Rmd` as `mill_v_err_abs` |
| `code/pdam_variance.R` | Calculates repeatability and measurement-error estimates for the *P. damicornis* case study; the resulting Blomqvist error ratio is carried into `RTM_analysis.Rmd` as `pdam_k` |
| `code/RTM_functions.R` | Stand-alone functions for the Blomqvist and no-heterogeneity method-of-moments estimators described in the Supplementary Text, provided to help users easily apply these methods to their own two-timepoint baseline-change datasets; included as companion utility code rather than as a direct step in the main analysis pipeline |

### Data

| File | Description |
|---|---|
| `data/PdamRbleaching.csv` | Case study data for the *P. damicornis* relationship between initial symbiont density and bleaching, sourced from [Cunning & Baker 2013](https://doi.org/10.1038/nclimate1711). Variables: `colony`, coral colony identifier; `sym`, dominant algal symbiont clade; `juntotal`, total symbiont density before bleaching, measured as symbiont:host cell ratio; `augtotal`, total symbiont density after bleaching, measured as symbiont:host cell ratio. |
| `data/PdamRwarming.csv` | Additional *P. damicornis* data used to estimate repeatability for the bleaching case study, sourced from [Cunning & Baker 2013](https://doi.org/10.1038/nclimate1711). Variables: `colony`, coral colony identifier; `sym`, dominant algal symbiont clade; `date`, sampling month; `total`, total symbiont density, measured as symbiont:host cell ratio. |
| `data/subMil1.csv` | Case study data for the *A. cervicornis* relationship between initial size and growth, derived from [Million et al. 2022](https://doi.org/10.1073/pnas.2203925119); associated data/code archive: [Zenodo](https://doi.org/10.5281/zenodo.7181336). Variables: `InitialSize`, initial total linear extension for the growth interval, in cm; `size`, final total linear extension for the growth interval, in cm; `diff`, growth interval label. |
| `data/Acer3DMorphologyData.csv` | Additional *A. cervicornis* data used to estimate measurement error for the growth case study, derived from [Million et al. 2022](https://doi.org/10.1073/pnas.2203925119); associated data/code archive: [Zenodo](https://doi.org/10.5281/zenodo.7181336). Variables: `Genotype`, coral genotype identifier; `FragID`, fragment identifier within genotype; `T0_TLE` and `T0_FieldTLE`, paired baseline measurements of total linear extension for the same fragment, in cm. |

### Figures

| File | Description |
|---|---|
| `figures/fig1.png` | Simulation examples showing how true slope, change heterogeneity, and measurement error affect observed and corrected baseline-change relationships |
| `figures/fig2.png` | Bias map comparing slope estimators across measurement-error levels and change heterogeneity |
| `figures/fig3.png` | Empirical tests of regression-to-the-mean corrections in the two coral case studies |

## Software

Analyses were run in R 4.3.1 with the following package versions: `tidyverse` 2.0.0, `ggpubr` 0.6.0, `ggpmisc` 0.6.1, `broom` 1.0.4, `patchwork` 1.3.0, `boot` 1.3-28.1, `future` 1.33.0, `furrr` 0.3.1, `ggtext` 0.1.2, `rptR` 0.9.23, `rmarkdown` 2.25, and `knitr` 1.45.
