
###############
## load data ##
###############
setwd("/Users/shiyaoxu/Documents/Research/causal_ODHSI/major_cardio_event")
d_treatment <- read.table("treatment.txt", header = FALSE)$V1
d_censoring_time <- read.table("censoring_time.txt", header = FALSE)$V1
d_event_time <- read.table("event_time.txt", header = FALSE)$V1
d_outcome <- read.table("outcome.txt", header = FALSE)$V1
d_outcome <- ifelse(d_outcome %in% c(2, 3, 4), 1, d_outcome)
d_time <- d_event_time
d_time[which(d_outcome == 0)] <- d_censoring_time[which(d_outcome == 0)]
rm(list=c("d_censoring_time", "d_event_time"))
covariates <- data.table::fread("sparse_design_matrix.txt")
covariates[which(covariates$j == 1662 & covariates$val == 1),]$i[90151] <- "1047"
covariates$i <- as.integer(covariates$i)
cov <- Matrix::sparseMatrix(i = covariates$i, j = covariates$j, x = covariates$val, repr = "T")

## parameters
rowId <- 1:length(d_outcome)
n <- length(id)



###################################
## PS calculation from real data ##
###################################

## calculate propensity score (no cross-fit)
treatProb <- estimateTreatProb(id=id, treatment=d_treatment, covariates=covariates, covIdTreatProb=NULL, 
                               treatProbEstimate="LASSO", maxCohortSizeForFitting=25000,
                               index_ls=NULL, crossFitNum=1)



###################################################
## survival and censoring hazards from real data ##
###################################################

## survival hazards and selected covariates
survHaz_initial <- estimateSimulationParams(treatment=d_treatment, covariates=covariates, 
                                            outcome=d_outcome, nInt=50, hazEstimate="survival", 
                                            hazMethod="ns", seed=2019)

## censoring hazards and selected covariates
cenHaz_initial <- estimateSimulationParams(treatment=d_treatment, covariates=covariates, 
                                            outcome=d_outcome, nInt=50, hazEstimate="censoring", 
                                            hazMethod="ns", seed=2019)


################
## simulation ##
################






## simulation
set.seed(2022)
SimData <- simData(survHaz=SurvHaz, cenHaz=CenHaz, treatment=d_treatment, 
                   maxTimeSurv=maxTimeSurv, maxTimeCen=maxTimeCen)

## counterfactuals
Truth <- counterFactuals(survHaz=SurvHaz, maxTime=maxTimeSurv)


## save results





########################################################
## survival and censoring hazards from simulated data ##
########################################################





########################
## estimate S(1)-S(0) ##
########################

## predicted hazards


## algorithm


## predicted hazards

## algorithm


####################
## plot S(1)-S(0) ##
####################






##############################
## estimate RMST(1)-RMST(0) ##
##############################





##########################
## plot RMST(1)-RMST(0) ##
##########################














