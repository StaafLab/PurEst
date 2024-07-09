<font size="12">**PurEst** </font>  <a href="https://staaflab.github.io/PurEst/"><img src="man/figures/logo.svg" align="right" height="278" alt="PurEst website" /></a>

Method for calculating tumor purity and adjusting Illumina 450/850K DNA methylation beta values. This work has been described in the article "**Tumor purity estimated from bulk DNA methylation can be used for adjusting beta values of individual samples to better reflect tumor biology**" published by [Sasiain, Nacer et al. in *JOURNAL*](LINK).

## How to install

```R
# Install devtools if required
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", quietly = TRUE) }

# Load the library
library(devtools)

# Install and load PurEst
install_github("StaafLab/PurEst/PurEst")
library(PurEst)
```

## Framework
The PurEst package is divided into four main functions:

1. *Correction of beta values from a cohort with known sample purities*: The function ***beta_correction_for_cohorts()*** performs beta value adjustment based on tumour sample compusition based on the original Staaf & Aine beta correction approch, correcting betas from a cohort of samples with known sample purities. More information about this approach can be found in the original publication: [Staaf & Aine, PLosOne, 2022](https://doi.org/10.1371/journal.pone.0265557).

2. *Creation of reference data from a cohort*: in the ***reference_regressions_generator()*** function, reference linear regressions are calculated based on DNA methylation beta values and tumor purity estimates of a cohort. Each regression represents a sample population as shown in the figure below. Scripts used in this module: ref_regression_calculator.r (main), new_function_correctBetas.r.
<div align="center">
  <img src="./man/figures/module1.png" width="350" class="center">
</div>
</p>

3. *Estimation of tumor purities for individual samples*: in thE ***PurEst()*** function, CpGs are filtered based on beta variance and then each CpG is processed individually per sample as shown in the figure below. Scripts used in this module: purity_estimator.r (main), predicting_purity.r, purity_coverage.r.

![](./man/figures/module2.png)

4. *Adjustment of beta values per CpG per sample*: in the ***reference_based_beta_correction()*** function, beta values are adjusted for tumor cells and inferred for normal cells using the reference regressions and the estimated purities. This can be carried out after refitting the regressions to include the new data points with estimated purities or using the original reference regressions. Scripts used in this module: final_beta_correction.r (main), new_function_correctBetas.r.

![](./man/figures/module3.png)
