#include <RcppEigen.h>
#include "CausalSurvival_types.h"
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd computeSubsetInformationMatrix(
  const Map<CsrMatrix> X, const Map<VectorXd> weight, int subsetSize
) {
  int nPred = X.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred);
  DiagonalMatrix subsetWeightMat = weight.head(subsetSize).asDiagonal();
  CsrMatrix subsetX = X.topRows(subsetSize);
  subsetInfoMat = subsetX.transpose() * subsetWeightMat * subsetX;
  return subsetInfoMat;
}



