# RegToMean

A simulation and analysis project for evaluating regression-to-the-mean (RTM) artifacts in relationships between trait changes and initial measurements, and methods to correct for them.


## 📂 Repository File Manifest  

### Main
| File | Description |
|---|---|
| `README.md` | Repository documentation and project overview |
| `RegToMean.Rproj` | RStudio project file for workspace management |

### Code
| File | Purpose |
|---|---|
| `code/RTM_analysis.Rmd` | Core analysis notebook — simulations, RTM correction method evaluation, figure generation |
| `code/RegressionToMean_forGithub.R` | Script for calculating measurement error for *A. cervicornis* case study |
| `code/pdam_variance.R` | Script for calculating measurement error for *P. damicornis* case study |

### Data
| File | Description |
|---|---|
| `data/PdamRbleaching.csv` | Case study data for *P. damicornis* releationship between initial symbiont density and bleaching |
| `data/PdamRwarming.csv` | Additional *P. damicornis* symbiont density timeseries data to compute repeatability |
| `data/subMil1.csv` | Case study data for *A. cervicornis* relationship between initial size and growth |

### Figures (generated outputs)
| File | Description |
|---|---|
| `figures/fig1.png` | Visualization of four canonical RTM scenarios comparing observed and corrected Δ vs. initial slopes under combinations of true slope, individual heterogeneity, and measurement error. |
| `figures/fig2.png` | Bias map showing how observed, Kelly–Price, and Blomqvist slope estimators deviate from the true slope across measurement-error levels and change heterogeneity. |
| `figures/fig3.png` | Empirical tests of RTM using coral bleaching + growth datasets showing how slope corrections change inference relative to naive regression. |

---
