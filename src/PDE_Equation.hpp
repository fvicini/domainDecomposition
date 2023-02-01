#ifndef __PDE_Equation_H
#define __PDE_Equation_H

#include "Eigen/Eigen"

namespace DOMAIN_DECOMPOSITION
{
  /// \brief 2D Primal Conforming Constant Lagrange Element Degree variable
  class PDE_Equation final
  {
    public:
      static Eigen::MatrixXd ComputeStiffnessMatrix(const unsigned int& numCellDofs,
                                                    const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                    const Eigen::VectorXd& quadratureWeights);

      static Eigen::VectorXd ComputeCellForcingTerm(const Eigen::VectorXd& forcingTermValues,
                                                    const Eigen::MatrixXd& basisFunctionValues,
                                                    const Eigen::VectorXd& quadratureWeights);

      static Eigen::VectorXd ForcingTerm(const Eigen::MatrixXd& points);
      static Eigen::VectorXd ExactSolution(const Eigen::MatrixXd& points);
      static Eigen::VectorXd ExactDerivativeSolution(const unsigned int& direction,
                                                     const Eigen::MatrixXd& points);
  };
}

#endif // __PDE_Equation_H
