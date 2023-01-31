#include "Fem2DSquareLagrangePCC.hpp"

using namespace std;
using namespace Eigen;

namespace DOMAIN_DECOMPOSITION
{
  // ***************************************************************************
  SquareMapping::Map SquareMapping::Compute(const Eigen::MatrixXd& vertices,
                                            const double& area) const
  {
    Map map;

    map.Q = Eigen::Matrix3d::Identity();
    map.Q(0, 0) = sqrt(area);
    map.Q(1, 1) = map.Q(0, 0);

    map.QInv = Eigen::Matrix3d::Identity();
    map.QInv(0, 0) = 1.0 / map.Q(0, 0);
    map.QInv(1, 1) = map.QInv(0, 0);

    map.detQ = area;
    map.b = vertices.col(0);

    return map;
  }
  // ***************************************************************************
  Eigen::MatrixXd Fem2DSquareLagrangePCC::EvaluateLambda(const MatrixXd& points) const
  {
    MatrixXd lambda = MatrixXd::Zero(points.cols(), 4);

    lambda.col(0) = 1.0 - points.row(0).array();
    lambda.col(1) = 1.0 - points.row(1).array();
    lambda.col(2) = points.row(0);
    lambda.col(3) = points.row(1);

    return lambda;
  }
  // ***************************************************************************
  vector<MatrixXd> Fem2DSquareLagrangePCC::EvaluateGradLambda(const MatrixXd& points) const
  {
    vector<MatrixXd> gradLambda(2, MatrixXd::Zero(points.cols(), 4));

    gradLambda[0].col(0).setConstant(-1.0);
    gradLambda[1].col(0).setZero();

    gradLambda[0].col(1).setZero();
    gradLambda[1].col(1).setConstant(-1.0);

    gradLambda[0].col(2).setOnes();
    gradLambda[1].col(2).setZero();

    gradLambda[0].col(3).setZero();
    gradLambda[1].col(3).setOnes();

    return gradLambda;
  }
  // ***************************************************************************
  Fem2DSquareLagrangePCC::LocalSpace Fem2DSquareLagrangePCC::Compute()
  {
    LocalSpace localSpace;

    localSpace.NumberBasisFunctions = 4;

    localSpace.ReferenceElementDofPositions.setZero(3, localSpace.NumberBasisFunctions);
    localSpace.ReferenceElementDofPositions.col(0)<< 0.0, 0.0, 0.0;
    localSpace.ReferenceElementDofPositions.col(1)<< 1.0, 0.0, 0.0;
    localSpace.ReferenceElementDofPositions.col(2)<< 1.0, 1.0, 0.0;
    localSpace.ReferenceElementDofPositions.col(3)<< 0.0, 1.0, 0.0;

    return localSpace;
  }
  // ***************************************************************************
  Eigen::MatrixXd Fem2DSquareLagrangePCC::Reference_BasisFunctions(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                                   const MatrixXd& points) const
  {
    MatrixXd values = MatrixXd::Zero(points.cols(),
                                     localSpace.NumberBasisFunctions);

    const MatrixXd lambda = EvaluateLambda(points);

    values.col(0) = lambda.col(0).array() * lambda.col(1).array();
    values.col(1) = lambda.col(2).array() * lambda.col(1).array();
    values.col(2) = lambda.col(2).array() * lambda.col(3).array();
    values.col(3) = lambda.col(0).array() * lambda.col(3).array();

    return values;
  }
  // ***************************************************************************
  vector<MatrixXd> Fem2DSquareLagrangePCC::Reference_BasisFunctionDerivatives(const Fem2DSquareLagrangePCC::LocalSpace& localSpace,
                                                                              const MatrixXd& points) const
  {
    vector<MatrixXd> values(2, MatrixXd::Zero(points.cols(),
                                              localSpace.NumberBasisFunctions));

    const MatrixXd lambda = EvaluateLambda(points);

    const vector<MatrixXd> gradLambda = EvaluateGradLambda(points);

    values[0].col(0) =
        gradLambda[0].col(0).array() * lambda.col(1).array() +
        lambda.col(0).array() * gradLambda[0].col(1).array();
    values[1].col(0) =
        gradLambda[1].col(0).array() * lambda.col(1).array() +
        lambda.col(0).array() * gradLambda[1].col(1).array();

    values[0].col(1) =
        gradLambda[0].col(2).array() * lambda.col(1).array() +
        lambda.col(2).array() * gradLambda[0].col(1).array();
    values[1].col(1) =
        gradLambda[1].col(2).array() * lambda.col(1).array() +
        lambda.col(2).array() * gradLambda[1].col(1).array();

    values[0].col(2) =
        gradLambda[0].col(2).array() * lambda.col(3).array() +
        lambda.col(2).array() * gradLambda[0].col(3).array();
    values[1].col(2) =
        gradLambda[1].col(2).array() * lambda.col(3).array() +
        lambda.col(2).array() * gradLambda[1].col(3).array();

    values[0].col(3) =
        gradLambda[0].col(0).array() * lambda.col(3).array() +
        lambda.col(0).array() * gradLambda[0].col(3).array();
    values[1].col(3) =
        gradLambda[1].col(0).array() * lambda.col(3).array() +
        lambda.col(0).array() * gradLambda[1].col(3).array();

    return values;
  }
  // ***************************************************************************
  std::vector<MatrixXd> Fem2DSquareLagrangePCC::Map_BasisFunctionDerivatives(const LocalSpace& localSpace,
                                                                             const SquareMapping::Map& map,
                                                                             const std::vector<Eigen::MatrixXd>& reference_values) const
  {
    std::vector<Eigen::MatrixXd> mapped_values(2,
                                               Eigen::MatrixXd::Zero(reference_values[0].rows(),
                                               localSpace.NumberBasisFunctions));

    for(unsigned int i = 0; i < 2; i++)
    {
      mapped_values[i] = map.QInv(i,i) *
                         reference_values[i];
      for(unsigned int j = 0; j < i; j++)
      {
        mapped_values[i] += map.QInv(j,i) *
                            reference_values[j];
        mapped_values[j] += map.QInv(i,j) *
                            reference_values[i];
      }
    }

    return mapped_values;
  }
  // ***************************************************************************
}
