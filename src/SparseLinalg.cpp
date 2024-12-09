#include <RcppEigen.h>
#include "CausalSurvival_types.h"
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd computeSubsetSparseInformationMatrix(
  const Map<CsrMatrix> X, const Map<VectorXd> weight, int subsetSize
) {
  int nPred = X.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred);
  DiagonalMatrix subsetWeightMat = weight.head(subsetSize).asDiagonal();
  CsrMatrix subsetX = X.topRows(subsetSize);
  subsetInfoMat = subsetX.transpose() * subsetWeightMat * subsetX;
  return subsetInfoMat;
}

// [[Rcpp::export]]
MatrixXd computeSubsetInformationMatrix(
    const Map<MatrixXd> Y, const Map<VectorXd> weight, int subsetSize
) {
  int nPred = Y.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred);
  DiagonalMatrix subsetWeightMat = weight.head(subsetSize).asDiagonal();
  MatrixXd subsetY = Y.topRows(subsetSize);
  subsetInfoMat = subsetY.transpose() * subsetWeightMat * subsetY;
  return subsetInfoMat;
}

// [[Rcpp::export]]
MatrixXd computeSubsetMatrix(
    const Map<CsrMatrix> X, const Map<VectorXd> weight, const Map<MatrixXd> Y, int subsetSize
) {
  int nPred = X.cols();
  int nPred2 = Y.cols();
  MatrixXd subsetMat = MatrixXd::Zero(nPred, nPred2);
  DiagonalMatrix subsetWeightMat = weight.head(subsetSize).asDiagonal();
  CsrMatrix subsetX = X.topRows(subsetSize);
  MatrixXd subsetY= Y.topRows(subsetSize);
  subsetMat = subsetX.transpose() * subsetWeightMat * subsetY;
  return subsetMat;
}

// [[Rcpp::export]]
VectorXd computeSubsetSparseMatVec(
    const Map<CsrMatrix> X, const Map<VectorXd> v, int subsetSize, bool transposed = true
) {
  CsrMatrix subsetX = X.topRows(subsetSize);
  VectorXd subsetv = v.head(subsetSize);
  if (transposed) {
    return subsetX.transpose() * subsetv;
  }else {
  return subsetX * v;
  }
}


// [[Rcpp::export]]
VectorXd computeSubsetMatVec(
    const Map<MatrixXd> Y, const Map<VectorXd> v, int subsetSize, bool transposed = true
) {
  MatrixXd subsetY = Y.topRows(subsetSize);
  VectorXd subsetv = v.head(subsetSize);
  if (transposed) {
    return subsetY.transpose() * subsetv;
  }else {
    return subsetY * v;
  }
}



