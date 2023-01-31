#ifndef __FEM2DSQUARELAGRANGEPCC_H
#define __FEM2DSQUARELAGRANGEPCC_H

#include "Eigen/Eigen"

namespace DOMAIN_DECOMPOSITION
{
  /// \brief Map square from reference [0,1]x[0,1] to generic square
  class SquareMapping
  {
    public:
      struct Map
      {
          Eigen::Matrix3d Q;
          Eigen::Matrix3d QInv;
          double detQ;
          Eigen::VectorXd b;
      };

    public:
      Map Compute(const Eigen::MatrixXd& vertices,
                  const double& area) const;

      inline Eigen::MatrixXd F(const Map& map,
                               const Eigen::MatrixXd& points) const
      { return (map.Q * points).colwise() + map.b; }
      inline Eigen::MatrixXd FInv(const Map& map,
                                  const Eigen::MatrixXd& points) const
      { return map.QInv * (points.colwise() - map.b); }

  };

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
      LocalSpace Compute();

      Eigen::MatrixXd Reference_BasisFunctions(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                               const Eigen::MatrixXd& points) const;

      std::vector<Eigen::MatrixXd> Reference_BasisFunctionDerivatives(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                                      const Eigen::MatrixXd& points) const;

      inline Eigen::MatrixXd Map_BasisFunctions(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                const SquareMapping::Map& map,
                                                const Eigen::MatrixXd& reference_values) const
      { return reference_values; }

      std::vector<Eigen::MatrixXd> Map_BasisFunctionDerivatives(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                                const SquareMapping::Map& map,
                                                                const std::vector<Eigen::MatrixXd>& reference_values) const;
  };
}

#endif // __FEM2DSQUARELAGRANGEPCC_H
