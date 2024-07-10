#' Estimate cross-fitted nuisance parameter: censoring hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param gR.estimate Model for estimating nuisance parameter: censoring hazards
#' @export

estimateNuisanceGR <- function(dlong, J, gR.estimate="glm"){

  ## container
  ID <- Time <- gR1 <- gR0 <- c()


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)


    ## model: glm
    fitR <- glm(Rt ~ treatment * (t + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = Jt == 1, family = binomial())


    ## predict
    gR1temp <- bound01(predict(fitR, newdata = mutate(d_test, treatment = 1), type = 'response'))
    gR0temp <- bound01(predict(fitR, newdata = mutate(d_test, treatment = 0), type = 'response'))


    ## store
    ID <- c(ID, idx_test)
    Time <- c(Time, d_test$t)
    gR1 <- c(gR1, gR1temp)
    gR0 <- c(gR0, gR0temp)
  }
  ## result
  out <- data.frame(id=ID, t=Time, gR1=gR1, gR0=gR0)
  return(out)
}


