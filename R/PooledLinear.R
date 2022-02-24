




coef_lm <- function(dlong, timeIntMidPoint, n, printIter){

    ## maxTime
    maxTime <- length(timeIntMidPoint)

    ## dlong
    dlong <- dlong[order(dlong$time, 1-dlong$eventObserved, decreasing = TRUE), ]

    ## subset index for each time point
    indx_subset <- sapply(timeIntMidPoint, function(x){sum(dlong$time[dlong$t == timeIntMidPoint[1]] >= x)}, USE.NAMES = FALSE)

    X_baseline <- cbind(rep(1, n), as(matrix(0, nrow = n, ncol = maxTime-1), "sparseMatrix"),
                        as(matrix(0, nrow = n, ncol = maxTime-1), "sparseMatrix"),
                        dlong$treat[dlong$t == timeIntMidPoint[1]])

    ## matrix
    Xy <- rep(0, length=dim(X_baseline)[2])
    XX <- matrix(0, nrow=dim(X_baseline)[2], ncol=dim(X_baseline)[2])

    for (i in 1:maxTime){

      if(i == 2){
        X_baseline <- cbind(rep(1, n), rep(1, n), as(matrix(0, nrow = n, ncol = maxTime-2), "sparseMatrix"),
                            dlong$treat[dlong$t == timeIntMidPoint[1]],
                            as(matrix(0, nrow = n, ncol = maxTime-2), "sparseMatrix"),
                            dlong$treat[dlong$t == timeIntMidPoint[1]])
      }else if(i > 2){
        X_baseline <- cbind(rep(1, n), as(matrix(0, nrow = n, ncol = i-2), "sparseMatrix"),
                            rep(1, n), as(matrix(0, nrow = n, ncol = maxTime-i), "sparseMatrix"),
                            as(matrix(0, nrow = n, ncol = i-2), "sparseMatrix"),
                            dlong$treat[dlong$t == timeIntMidPoint[1]],
                            as(matrix(0, nrow = n, ncol = maxTime-i), "sparseMatrix"),
                            dlong$treat[dlong$t == timeIntMidPoint[1]])
      }

      X_baseline <- Matrix::sparseMatrix(i = Matrix::summary(X_baseline)$i, j = Matrix::summary(X_baseline)$j,
                                         x = Matrix::summary(X_baseline)$x, repr = "R")

      ## Xy
      Xy_temp <- computeSubsetSparseMatVec(X=X_baseline, v=dlong$Lt[which(dlong$t == timeIntMidPoint[i])],
                                           subsetSize=indx_subset[i], transposed=TRUE)
      Xy <- Xy + Xy_temp
      ## XX
      XX_temp <- computeSubsetSparseInformationMatrix(X=X_baseline, weight=rep(1, n), subsetSize=indx_subset[i])
      XX <- XX + XX_temp

      if(printIter){print(i)}

    }

    ## beta
    beta <- solve(XX, Xy)

    ## predict
    X_baseline_predict <- expand.grid(t=1:maxTime, treat=1)
    X_baseline_predict <- as.matrix(cbind(model.matrix(~ as.factor(t), X_baseline_predict), model.matrix(~ as.factor(t), X_baseline_predict)[, -1], X_baseline_predict[,"treat"]))
    Haz1 <- X_baseline_predict %*% beta

    X_baseline_predict <- expand.grid(t=1:maxTime, treat=0)
    X_baseline_predict <- as.matrix(cbind(model.matrix(~ as.factor(t), X_baseline_predict), model.matrix(~ as.factor(t), X_baseline_predict)[, -1]*0, X_baseline_predict[,"treat"]))
    Haz0 <- X_baseline_predict %*% beta

    ## result
    out <- data.frame(Haz1=Haz1, Haz0=Haz0)
    return(out)
}



