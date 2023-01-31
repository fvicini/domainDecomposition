#ifndef __FEM2DSQUARELAGRANGEPCC_H
#define __FEM2DSQUARELAGRANGEPCC_H

#include "Eigen/Eigen"

namespace GeDiM
{
  /// \brief 2D Primal Conforming Constant Lagrange Element Degree variable
  class Fem2DSquareLagrangePCC final
  {
    public:
      struct LocalSpace
      {
          unsigned int NumberBasisFunctions; ///< Number of total basis functions
          Eigen::MatrixXd ReferenceElementDofPositions; ///< reference element dof points
      };

    private:
      /// \brief evaluate the lambda function for basis functions (Barycentric coordinate system)
      Eigen::MatrixXd EvaluateLambda(const Eigen::MatrixXd& points) const;
      /// \brief evaluate the lambda gradient function for basis functions
      std::vector<Eigen::MatrixXd> EvaluateGradLambda(const Eigen::MatrixXd& points) const;

    public:
      Fem2DSquareLagrangePCC();
      ~Fem2DSquareLagrangePCC();

      LocalSpace Initialize();

      Eigen::MatrixXd EvaluateBasisFunctions(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                             const Eigen::MatrixXd& points) const;

      std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                                    const Eigen::MatrixXd& points) const;
  };
}

#endif // __FEM2DSQUARELAGRANGEPCC_H
