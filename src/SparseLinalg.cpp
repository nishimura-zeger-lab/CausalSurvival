#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map; // 'Map' (i.e. reference without making copies) R matrices
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> EigenDiagonalMatrix;
typedef Map<SparseMatrix<double, Eigen::RowMajor>> MappedCsr;


// [[Rcpp::export]]
MatrixXd computeSubsetInformationMatrix(
  const MappedCsr X, const VectorXd sqrtWeight, int subsetSize
) {
  int nPred = X.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred);\
  EigenDiagonalMatrix subsetSqrtWeightMat = sqrtWeight.head(subsetSize).asDiagonal();
  SparseMatrix<double, Eigen::RowMajor> subsetWeigtedX = subsetSqrtWeightMat * X.topRows(subsetSize);
  subsetInfoMat = subsetWeigtedX.transpose() * subsetWeigtedX;
  return subsetInfoMat;
}



