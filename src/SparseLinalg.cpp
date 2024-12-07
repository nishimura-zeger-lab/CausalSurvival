#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map; // 'Map' (i.e. reference without making copies) R matrices
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Map<SparseMatrix<double, Eigen::ColMajor>> MappedCsc;
typedef Map<SparseMatrix<double, Eigen::RowMajor>> MappedCsr;

// [[Rcpp::export]]
MatrixXd computeSubsetInformationMatrix(
  const MappedCsr X, const Map<VectorXd> weight, int subsetSize
) {
  int nPred = X.cols();
  MatrixXd subsetInfoMat = MatrixXd::Zero(nPred, nPred);
  subsetInfoMat = X.topRows(subsetSize).transpose() * weight.head(subsetSize).asDiagonal() * X.topRows(subsetSize);
  return subsetInfoMat;
}



