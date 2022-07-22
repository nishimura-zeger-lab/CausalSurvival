# CausalSurvival: causal survival methods for large-scale observational data

The CausalSurvival R package implements several causal survival methods for analyzing large-scale observational data with point treatment, baseline covariates and discrete time-to-event outcome with informative right-censoring. The estimation methods include inverse probability weighting (IPW, weighting by the inverse of propensity score and censoring probability), propensity score stratification (stratified Cox model), Targeted Maximum Likelihood Estimation (TMLE), and augmented IPW (AIPW) in estimating the difference in survival probability and restricted mean survival time at one or more pre-specified time points. 

### Install package from GitHub

``` r
remotes::install_github("nishimura-zeger-lab/CausalSurvival")
```

## Estimation methods

This section shows how to run the algorithms. First, we need to estimate the nuisance parameters (propensity score, survival and censoring hazards) with LASSO logistic regression. 

### Estimate nuisance parameters

```r

## propensity score
ps <- estimateTreatProb(id=id, treatment=treatment, covariates=covariates)

## survival hazards
estimatedSurvHazard <- estimateHazards(coarsenedData=coarsenedData, outcome=outcome, 
                                       treatment=treatment, covariates=covariates, 
                                       covId=survCov_indx,
                                       hazEstimate="survival", hazMethod="ns")
                                       
## censoring hazards
estimatedCenHazard <- estimateHazards(coarsenedData=coarsenedData, outcome=outcome, 
                                      treatment=treatment, covariates=covariates, 
                                      covId=cenCov_indx,
                                      hazEstimate="censoring", hazMethod="ns")
                                      
```

### IPW

```r

ipw <- computeEstimator(treatment=treatment, outcome=outcome, time=coarsenedData$timeInt,
                       initial_survHaz=NULL, cenHaz=estimatedCenHazard$haz, treatProb=treatProb,
                       coarseningParam=coarsenedData, method="IPW")

```

### Stratified Cox model

```r

cox <- CausalSurvival::computeEstimator(treatment=treatment, outcome=outcome, time=coarsenedData$timeInt,
                                        initial_survHaz=estimatedSurvHazard$haz, cenHaz=estimatedCenHazard$haz, 
                                        treatProb=treatProb, coarseningParam=coarsenedData_sim, 
                                        method="cox", hazMethod="ns")

weightedCox <- CausalSurvival::computeEstimator(treatment=treatment, outcome=outcome, time=coarsenedData$timeInt,
                                                initial_survHaz=estimatedSurvHazard$haz, cenHaz=estimatedCenHazard$haz, 
                                                treatProb=treatProb, coarseningParam=coarsenedData_sim, 
                                                method="weightedCox", hazMethod="ns")

```

### TMLE

```r

tmle <- computeEstimator(treatment=treatment, outcome=outcome, time=coarsenedData$timeInt,
                         initial_survHaz=estimatedSurvHazard$haz, cenHaz=estimatedCenHazard$haz, treatProb=treatProb,
                         coarseningParam=coarsenedData, method="TMLE")

```

### AIPW

```r

aipw <- computeEstimator(treatment=treatment, outcome=outcome, time=coarsenedData$timeInt,
                         initial_survHaz=estimatedSurvHazard$haz, cenHaz=estimatedCenHazard$haz, treatProb=treatProb,
                         coarseningParam=coarsenedData, method="AIPW")

```

More information on using estimation methods with a simulated dataset can be found in Vignette. 

## References










