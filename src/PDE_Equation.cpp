#include "PDE_Equation.hpp"

using namespace std;
using namespace Eigen;

namespace DOMAIN_DECOMPOSITION
{
  // ***************************************************************************
  Eigen::MatrixXd PDE_Equation::ComputeStiffnessMatrix(const unsigned int& numCellDofs,
                                                       const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    MatrixXd cellMatrixA;
    cellMatrixA.setZero(numCellDofs,
                        numCellDofs);

    for(unsigned int i = 0; i < basisFunctionDerivativeValues.size(); i++)
      cellMatrixA += basisFunctionDerivativeValues[i].transpose() *
                     quadratureWeights.asDiagonal() *
                     basisFunctionDerivativeValues[i];

    return cellMatrixA;
  }
  // ***************************************************************************
  Eigen::VectorXd PDE_Equation::ComputeCellForcingTerm(const Eigen::VectorXd& forcingTermValues,
                                                       const Eigen::MatrixXd& basisFunctionValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    return
        basisFunctionValues.transpose() *
        quadratureWeights.asDiagonal() *
        forcingTermValues;
  }
  // ***************************************************************************
  Eigen::VectorXd ForcingTerm(const Eigen::MatrixXd& points)
  {
    return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                   points.row(0).array() * (1.0 - points.row(0).array()));
  }
  // ***************************************************************************
}
