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
  const MappedCsr X, const Map<VectorXd> weight, int subsetSize
) {
  int nPred = X.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred).selfadjointView<Eigen::Lower>();
  EigenDiagonalMatrix subsetWeightMat = weight.head(subsetSize).asDiagonal();
  SparseMatrix<double, Eigen::RowMajor> subsetX = X.topRows(subsetSize);
  subsetInfoMat = subsetX.transpose() * subsetWeightMat * subsetX;
  return subsetInfoMat;
}



