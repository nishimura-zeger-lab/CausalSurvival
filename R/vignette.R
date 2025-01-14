######################################
## method chosen for the simulation ##
######################################

nInt <- 50
hazMethod <- "ns"


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
survHaz_initial <- estimateSimulationParams(outcome=d_outcome, time=d_time, treatment=d_treatment, covariates=covariates,
                                            simOutcome=NULL, simTime=NULL, cov_indx=NULL,
                                            nInt=nInt, hazEstimate="survival", hazMethod=hazMethod, seed=2019)

## censoring hazards and selected covariates
cenHaz_initial <- estimateSimulationParams(outcome=d_outcome, time=d_time, treatment=d_treatment, covariates=covariates,
                                           simOutcome=NULL, simTime=NULL, cov_indx=NULL,
                                           nInt=nInt, hazEstimate="censoring", hazMethod=hazMethod, seed=2019)


################
## simulation ##
################

## counterfactuals
Truth <- counterFactuals(time=d_time, outcome=d_outcome, survHaz=survHaz_initial$haz, nInt=nInt)

## simulation
set.seed(2019)
SimData <- simData(time=d_time, outcome=d_outcome, treatment=d_treatment,
                   survHaz=survHaz_initial$haz, cenHaz=cenHaz_initial$haz, nInt=nInt)


########################################################
## survival and censoring hazards from simulated data ##
########################################################


## survival hazards and selected covariates
survHaz_initial <- estimateSimulationParams(outcome=d_outcome, time=d_time, treatment=d_treatment, covariates=covariates,
                                            simOutcome=SimData$ObservedEvent, simTime=SimData$ObservedTime, covId=survHaz_initial$cov_indx,
                                            nInt=nInt, hazEstimate="survival", hazMethod=hazMethod, seed=2019)

## censoring hazards and selected covariates
cenHaz_initial <- estimateSimulationParams(outcome=d_outcome, time=d_time, treatment=d_treatment, covariates=covariates,
                                           simOutcome=SimData$ObservedEvent, simTime=SimData$ObservedTime, covId=cenHaz_initial$cov_indx,
                                           nInt=nInt, hazEstimate="censoring", hazMethod=hazMethod, seed=2019)



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














