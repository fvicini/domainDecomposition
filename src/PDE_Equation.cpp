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
  VectorXd PDE_Equation::ForcingTerm(const Eigen::MatrixXd& points)
  {
    return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                   points.row(0).array() * (1.0 - points.row(0).array()));
  }
  // ***************************************************************************
  Eigen::VectorXd PDE_Equation::ExactSolution(const Eigen::MatrixXd& points)
  {
    return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                   points.row(0).array() * (1.0 - points.row(0).array()));
  }
  // ***************************************************************************
  Eigen::VectorXd PDE_Equation::ExactDerivativeSolution(const unsigned int& direction,
                                                        const Eigen::MatrixXd& points)
  {
    if (direction == 0)
      return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
    else if (direction == 1)
      return 16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array());
    else if (direction == 2)
      return Eigen::VectorXd::Zero(points.cols());
    else
      throw std::runtime_error("Error on direction");
  }
  // ***************************************************************************
}
