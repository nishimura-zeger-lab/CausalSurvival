#' Simulate time-to-event data with survival and censoring hazards
#'
#' @param survHaz Output from estimateHaz
#' @param cenHaz Output from estimateHaz
#' @param maxTimeSurv maximum time for event
#' @param maxTimeCen maximum time for censoring
#'

simData <- function(survHaz, cenHaz, treatment, maxTimeSurv, maxTimeCen){

  survHaz_all <- survHaz$Haz1 * treatment[survHaz$ID] + survHaz$Haz0 * (1 - treatment[survHaz$ID])
  cenHaz_all <- cenHaz$Haz1 * treatment[cenHaz$ID] + cenHaz$Haz0 * (1 - treatment[cenHaz$ID])

  rS <- rbinom(length(survHaz_all), 1, survHaz_all)
  rG <- rbinom(length(cenHaz_all), 1, cenHaz_all)

  tS <- tapply(rS, rep(1:n, each=maxTimeSurv), function(x){which(x == 1)[1]})
  tG <- tapply(rG, rep(1:n, each=maxTimeCen), function(x){which(x == 1)[1]})
  tG <- ifelse(is.na(tG), maxTimeCen, tG)

  ObservedTime <- apply(cbind(tS, tG), 1, function(x){min(x, na.rm = TRUE)})
  ObservedEvent <- 1*(tS<tG)
  ObservedEvent <- ifelse(is.na(ObservedEvent), 0, ObservedEvent)

  ObservedEvent <- as.double(ObservedEvent)
  ObservedTime <- as.double(ObservedTime)

  return(list(ObservedEvent=ObservedEvent,
              ObservedTime=ObservedTime))

}

#' Calculate counterfactuals from simulated data
#'
#'

counterFactuals <- function(survHaz, maxTime){

  ## paramters
  n <- dim(survHaz)[1]/maxTime

  ## surv prob
  Sm0 <- unlist(tapply(1-survHaz$Haz0, rep(1:n, each=maxTime), cumprod, simplify = FALSE), use.names = FALSE)
  S0 <- tapply(Sm0, rep(1:maxTime, n), mean)
  rm(list=c("Sm0"))
  Sm1 <- unlist(tapply(1-survHaz$Haz1, rep(1:n, each=maxTime), cumprod, simplify = FALSE), use.names = FALSE)
  S1 <- tapply(Sm1, rep(1:maxTime, n), mean)
  rm(list=c("Sm1"))

  ## rmst
  rmst0 <- cumsum(S0)
  rmst1 <- cumsum(S1)

  return(data.frame(S0=S0, S1=S1, rmst0=rmst0, rmst1=rmst1))

}








