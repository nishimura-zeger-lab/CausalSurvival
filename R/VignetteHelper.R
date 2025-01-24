#' Preprocess cohortMethodData
#'
#' @export
filterAndTidyCovariatesForPs <- function(cohortMethodData,
                                         population,
                                         excludeCovariateIds = c(),
                                         includeCovariateIds = c()) {

  covariates <- cohortMethodData$covariates %>%
    dplyr::filter(.data$rowId %in% local(population$rowId))
  if (length(includeCovariateIds) != 0) {
    covariates <- covariates %>%
      dplyr::filter(.data$covariateId %in% includeCovariateIds)
  }
  if (length(excludeCovariateIds) != 0) {
    covariates <- covariates %>%
      dplyr::filter(!.data$covariateId %in% excludeCovariateIds)
  }
  filteredCovariateData <- Andromeda::andromeda(covariates = covariates,
                                                covariateRef = cohortMethodData$covariateRef,
                                                analysisRef = cohortMethodData$analysisRef)
  metaData <- attr(cohortMethodData, "metaData")
  metaData$populationSize <- nrow(population)
  attr(filteredCovariateData, "metaData") <- metaData
  class(filteredCovariateData) <- "CovariateData"

  covariateData <- FeatureExtraction::tidyCovariateData(filteredCovariateData)
  Andromeda::close(filteredCovariateData)
  return(covariateData)
}
