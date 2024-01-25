**USAGE EXAMPLE**

This file explains the usage of the complete PurEst beta correction pipeline for samples without known purity values. The example dataset used consists of the 5000 most variable CpGs from the TCGA-BRCA data set, which was splitted into a subset containing 80% of the samples to generate the reference regressions, and the remaining 20% to estimate sample purity values and to correct betas.

* IMPORTANT: This example data set is not intended to be used to evaluate the tool's performace. Due to the low number of CpGs included aiming to speed up the process the performance may not be optimal.

The whole pipeline can be run following the next steps;

1. Generation of reference regressions
```bash
cd ./generating_regressions; 

Rscript ../../scripts/ref_regression_calculator.r -c 6 -b ../data/reference_data/betas_ref.rds -p ../data/reference_data/purity_ref.rds -o ./ -n example_ref;
```

2. Sample purity estimation

```bash
cd ../estimating_purities;

#Running the purity estimation script using 6 cores
Rscript ../../scripts/purity_estimator.r -c 6 -a 0.7 -s 0.25 -v 0.05 -p 0.96 -d ../generating_regressions/ -b ../data/data_to_correct/betas_toCorrect.rds -o example_estimated_purity;
```

3. Beta correction

```bash
#Correcting betas refitting the reference regressions
cd ../correcting_betas/refitting;

Rscript ../../../scripts/final_beta_correction.r -c 6 -F FALSE -o ./ -n example_refitting -r TRUE -P ../../data/reference_data/purity_ref.rds -B ../../data/reference_data/betas_ref.rds -b ../../data/data_to_correct/betas_toCorrect.rds -p ../../estimating_purities/example_estimated_purity.tsv;

#Correcting betas without refitting the reference regressions
cd ../not_refitting;

Rscript ../../../scripts/final_beta_correction.r -c 6 -F FALSE -o ./ -n example_refitting -r FALSE -R ../../generating_regressions -b ../../data/data_to_correct/betas_toCorrect.rds -p ../../estimating_purities/example_estimated_purity.tsv
``````
