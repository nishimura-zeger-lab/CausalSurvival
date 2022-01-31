#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map; // 'Map' (i.e. reference without making copies) R matrices
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrix;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> CsrMatrix;
typedef Map<SparseMatrix<double, Eigen::RowMajor>> MappedCsr;


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



