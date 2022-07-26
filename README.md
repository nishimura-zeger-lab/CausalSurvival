# CausalSurvival: causal survival methods for large-scale observational data

The CausalSurvival R package implements several causal survival methods for analyzing large-scale observational data with point treatment, baseline covariates and discrete time-to-event outcome with informative right-censoring. The estimation methods include inverse probability weighting (IPW, weighting by the inverse of propensity score and censoring probability), propensity score stratification (stratified Cox model), Targeted Maximum Likelihood Estimation (TMLE), and augmented IPW (AIPW) in estimating the difference in survival probability and restricted mean survival time at one or more pre-specified time points. 

### Install package from GitHub

``` r
remotes::install_github("nishimura-zeger-lab/CausalSurvival")
```

## Estimation methods

This section shows how to run the algorithms with the breast cancer dataset in the `survival` package. 

### Load dataset

To run the algorithms, we need to organize data into several data frames: id, treatment, outcome, time and covariates. We then further coarsen survival data into time intervals with each interval containing 10 events. 

```r
library(survival)
library(tidyr)
library(dplyr)
library(CausalSurvival)

data(cancer, package="survival")

id <- 1:dim(gbsg)[1]
treatment <- gbsg$hormon
outcome <- gbsg$status
time <- gbsg$rfstime

gbsg$size <- (gbsg$size-mean(gbsg$size))/sd(gbsg$size)
gbsg$pgr <- (gbsg$pgr-mean(gbsg$pgr))/sd(gbsg$pgr)
gbsg$er <- (gbsg$er-mean(gbsg$er))/sd(gbsg$er)
gbsg$grade <- ifelse(gbsg$grade == 1, 1, 0)
gbsg$nodes <- ifelse(gbsg$nodes == 1, 1, 0)
covariates <- gbsg %>% mutate(rowId = id) %>% 
                          gather(covariateId, covariateValue, -c(hormon, status, rfstime, rowId)) %>% 
                          select(-c(hormon, status, rfstime))
covariates$covariateId <- as.numeric(as.character(factor(covariates$covariateId, labels=1:8)))

coarsenedData <- coarsenData(time=time, outcome=outcome)

```

### Estimate nuisance parameters

First, estimate the nuisance parameters (propensity score, survival and censoring hazards) with ridge logistic regression. 

```r

## propensity score
estimatedPS <- estimateTreatProb(id=id, treatment=treatment, covariates=covariates)

## survival hazards
estimatedSurvHazard <- estimateHazards(coarsenedData=coarsenedData, outcome=outcome, 
                                       treatment=treatment, covariates=covariates, 
                                       hazEstimate="survival", hazMethod="ns")
                                       
## censoring hazards
estimatedCenHazard <- estimateHazards(coarsenedData=coarsenedData, outcome=outcome, 
                                      treatment=treatment, covariates=covariates, 
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

More information on using estimation methods with a simulated dataset to perform a simulation study and output the two types of causal estimands of interest specified in the protocol can be found in Vignette. 

## References










