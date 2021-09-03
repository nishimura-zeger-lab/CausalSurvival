#' Estimate cross-fitted nuisance parameter survival hazards
#'
#' @param dlong Long-format survival data from function transformData(dwide, freq.time)
#' @param J For cross-fitting
#' @param h.estimate Model for estimating nuisance parameter: survival hazards
#' @export

estimateNuisanceH <- function(dlong, J, h.estimate="glm"){

  ## container
  ID <- Time <- h1 <- h0 <- c()


  for (i in 1:J){
    ## training and testing sets
    idx_test <- index_ls[[i]]
    idx_train <- setdiff(unique(dlong$id), idx_test)
    d_test <- subset(dlong, id %in% idx_test)
    d_train <- subset(dlong, id %in% idx_train)


    ## Survival hazard: glm
    fitL <- glm(Lt ~ treatment * (t + age65 + cardiovascular + female + CHADS2),
                data = d_train, subset = It == 1, family = binomial())


    ## predict
    h1temp <- bound01(predict(fitL, newdata = mutate(d_test, treatment = 1), type = 'response'))
    h0temp <- bound01(predict(fitL, newdata = mutate(d_test, treatment = 0), type = 'response'))


    ## store
    ID <- c(ID, idx_test)
    Time <- c(Time, d_test$t)
    h1 <- c(h1, h1temp)
    h0 <- c(h0, h0temp)
  }
  ## result
  out <- data.frame(id=ID, t=Time, h1=h1, h0=h0)
  return(out)
}


