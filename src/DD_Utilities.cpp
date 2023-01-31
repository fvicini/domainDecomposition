#include "DD_Utilities.hpp"

#include "Fem2DSquareLagrangePCC.hpp"
#include "PDE_Equation.hpp"
#include "VTKUtilities.hpp"
#include "MapQuadrilateral.hpp"
#include "Quadrature_Gauss2D_Square.hpp"

#include "mpi.h"

using namespace std;

namespace DOMAIN_DECOMPOSITION
{
  // ***************************************************************************
  DD_Utilities::DD_Utilities()
  {
  }
  DD_Utilities::~DD_Utilities()
  {
  }
  // ***************************************************************************
  void DD_Utilities::PrintMessage(const int& rank,
                                  std::ostream& output,
                                  const std::string& message,
                                  const bool& onlyMaster)
  {
    if (onlyMaster && rank != 0)
      return;

    output<< ">> "<< "Process "<< rank<< ": "<< message<< endl;
  }
  // ***************************************************************************
  void DD_Utilities::Assert(const int& rank,
                            const bool& result,
                            const string& message)
  {
    if (!result)
    {
      if (!message.empty())
        PrintMessage(rank, cerr, message, false);

      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  // ***************************************************************************
  DD_Utilities::Point_Info DD_Utilities::Point_Info_Global(const unsigned int& i,
                                                           const unsigned int& j,
                                                           const unsigned int& n_1D_points)
  {
    Point_Info info;

    info.Index = n_1D_points * j + i;

    return info;
  }
  // ***************************************************************************
  DD_Utilities::Square_Info DD_Utilities::Square_Info_Global(const unsigned int& i,
                                                             const unsigned int& j,
                                                             const unsigned int& n_1D_points,
                                                             const unsigned int& n_1D_squares)
  {
    Square_Info info;

    info.Index = n_1D_squares * j + i;
    info.Points_Index.resize(4);
    info.Points_Index[0] = Point_Info_Global(i, j, n_1D_points).Index;
    info.Points_Index[1] = Point_Info_Global(i + 1, j, n_1D_points).Index;
    info.Points_Index[2] = Point_Info_Global(i + 1, j + 1, n_1D_points).Index;
    info.Points_Index[3] = Point_Info_Global(i, j + 1, n_1D_points).Index;

    return info;
  }

  // ***************************************************************************
  DD_Utilities::Domain_Info DD_Utilities::GetDomain_Info_Global(const int& rank,
                                                                const unsigned int& i_domain,
                                                                const unsigned int& j_domain,
                                                                const unsigned int& n_1D_domains,
                                                                const unsigned int& n_1D_squares_domain)
  {
    Domain_Info info;

    info.Axes_Index.resize(2);
    info.Axes_Index[0] = (rank % n_1D_domains) * n_1D_squares_domain + i_domain;
    info.Axes_Index[1] = (rank / n_1D_domains) * n_1D_squares_domain + j_domain;

    return info;
  }
  // ***************************************************************************
  DD_Utilities::Problem_Info DD_Utilities::ComputeProblemInfo(const int& rank,
                                                              const int& n_domains,
                                                              const Gedim::IMeshDAO& globalMesh)
  {
    Problem_Info info;
    info.Num_1D_domains = sqrt(n_domains);

    info.Num_2D_points = globalMesh.Cell0DTotalNumber();
    info.Num_2D_squares = globalMesh.Cell2DTotalNumber();
    info.Num_1D_points = sqrt(info.Num_2D_points);
    info.Num_1D_squares = info.Num_1D_points - 1;

    Assert(rank,
           info.Num_1D_squares % 2 == 0,
           "n_1D_squares shall be even");

    info.Num_1D_squares_domain = info.Num_1D_squares / info.Num_1D_domains;
    info.Num_1D_points_domain = info.Num_1D_squares_domain + 1;

    return info;
  }
  // ***************************************************************************
  void DD_Utilities::ExportDomainToVtu(const int& rank,
                                       const Problem_Info& problem_info,
                                       const Gedim::IMeshDAO& globalMesh,
                                       const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      for (unsigned int i_domain = 0; i_domain < problem_info.Num_1D_points_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < problem_info.Num_1D_points_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                problem_info.Num_1D_domains,
                                                                problem_info.Num_1D_squares_domain);

          const Point_Info point_info = Point_Info_Global(domain_info.Axes_Index.at(0),
                                                          domain_info.Axes_Index.at(1),
                                                          problem_info.Num_1D_points);

          const vector<double> cell0DDomain(1, rank);
          const vector<double> cell0DGlobalIndex(1, point_info.Index);
          exporter.AddPoint(globalMesh.Cell0DCoordinates(point_info.Index),
                            {
                              {
                                "Domain",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(cell0DDomain.size()),
                                cell0DDomain.data()
                              },
                              {
                                "GlobalIndex",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(cell0DGlobalIndex.size()),
                                cell0DGlobalIndex.data()
                              }
                            });
        }
      }

      exporter.Export(exportFolder +
                      "/Domain_" +
                      to_string(rank) +
                      "_Mesh_Cell0Ds" +
                      ".vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      for (unsigned int i_domain = 0; i_domain < problem_info.Num_1D_squares_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < problem_info.Num_1D_squares_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                problem_info.Num_1D_domains,
                                                                problem_info.Num_1D_squares_domain);

          const Square_Info square_info = Square_Info_Global(domain_info.Axes_Index.at(0),
                                                             domain_info.Axes_Index.at(1),
                                                             problem_info.Num_1D_points,
                                                             problem_info.Num_1D_squares);

          const vector<double> cell2DDomain(1, rank);
          const vector<double> cell2DGlobalIndex(1, square_info.Index);
          exporter.AddPolygon(globalMesh.Cell2DVerticesCoordinates(square_info.Index),
                              {
                                {
                                  "Domain",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(cell2DDomain.size()),
                                  cell2DDomain.data()
                                },
                                {
                                  "GlobalIndex",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(cell2DGlobalIndex.size()),
                                  cell2DGlobalIndex.data()
                                }
                              });
        }
      }

      exporter.Export(exportFolder +
                      "/Domain_" +
                      to_string(rank) +
                      "_Mesh_Cell2Ds" +
                      ".vtu");
    }
  }
  // ***************************************************************************
  DD_Utilities::DOF_Info DD_Utilities::CreateDOFs(const int& rank,
                                                  const int& n_domains,
                                                  const Problem_Info& problem_info,
                                                  const Gedim::IMeshDAO& globalMesh)
  {
    DOF_Info info;

    // initialize
    info.Num_Dirichlets = 0; // all external borders are dirichlets
    info.Num_Internals = 0;
    info.Num_Gamma = 0;
    info.Num_Globals = globalMesh.Cell0DTotalNumber();
    info.Cell0Ds_DOF.resize(globalMesh.Cell0DTotalNumber());

    // check mesh points type
    for (unsigned int i = 0; i < problem_info.Num_1D_points; i++)
    {
      for (unsigned int j = 0; j < problem_info.Num_1D_points; j++)
      {
        const Point_Info global_info = Point_Info_Global(i, j, problem_info.Num_1D_points);

        if (i == 0 || i == problem_info.Num_1D_points - 1 ||
            j == 0 || j == problem_info.Num_1D_points - 1)
        {
          info.Cell0Ds_DOF[global_info.Index].Type = DOF_Info::DOF::Types::Dirichlet;
          info.Num_Dirichlets++;
        }
        else if (i % problem_info.Num_1D_squares_domain == 0 ||
                 j % problem_info.Num_1D_squares_domain == 0)
        {
          info.Cell0Ds_DOF[global_info.Index].Type = DOF_Info::DOF::Types::Gamma;
          info.Num_Gamma++;
        }
        else
        {
          info.Cell0Ds_DOF[global_info.Index].Type = DOF_Info::DOF::Types::Internal;
          info.Num_Internals++;
        }
      }
    }

    // numerate dirichlet and gamma points
    unsigned int dirichlet_counter = 0;
    unsigned int gamma_counter = 0;
    for (unsigned int p = 0; p < globalMesh.Cell0DTotalNumber(); p++)
    {
      switch (info.Cell0Ds_DOF[p].Type)
      {
        case DOF_Info::DOF::Types::Dirichlet:
          info.Cell0Ds_DOF[p].GlobalIndex = info.Num_Internals +
                                            info.Num_Gamma +
                                            dirichlet_counter++;
          break;
        case DOF_Info::DOF::Types::Gamma:
          info.Cell0Ds_DOF[p].GlobalIndex = info.Num_Internals +
                                            gamma_counter++;
          break;
        case DOF_Info::DOF::Types::Internal:
          continue;
        default:
          throw runtime_error("unkwnon point type");
      }
    }

    Assert(rank,
           dirichlet_counter == info.Num_Dirichlets);
    Assert(rank,
           gamma_counter == info.Num_Gamma);

    // numerate internal dofs with processes order
    unsigned int internal_counter = 0;
    for (unsigned int rank = 0; rank < n_domains; rank++)
    {
      for (unsigned int j_domain = 0; j_domain < problem_info.Num_1D_points_domain; j_domain++)
      {
        for (unsigned int i_domain = 0; i_domain < problem_info.Num_1D_points_domain; i_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                problem_info.Num_1D_domains,
                                                                problem_info.Num_1D_squares_domain);

          const Point_Info point_info = Point_Info_Global(domain_info.Axes_Index.at(0),
                                                          domain_info.Axes_Index.at(1),
                                                          problem_info.Num_1D_points);

          switch (info.Cell0Ds_DOF[point_info.Index].Type)
          {
            case DOF_Info::DOF::Types::Dirichlet:
            case DOF_Info::DOF::Types::Gamma:
              continue;
            case DOF_Info::DOF::Types::Internal:
              info.Cell0Ds_DOF[point_info.Index].GlobalIndex = internal_counter++;
              break;
            default:
              throw runtime_error("unkwnon point type");
          }
        }
      }
    }

    Assert(rank,
           internal_counter == info.Num_Internals);

    return info;
  }
  // ***************************************************************************
  void DD_Utilities::ExportDOFsToVtu(const DOF_Info& dofs,
                                     const Gedim::IMeshDAO& globalMesh,
                                     const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      vector<double> dofsType(globalMesh.Cell0DTotalNumber());
      vector<double> dofsIndex(globalMesh.Cell0DTotalNumber());

      for (unsigned int p = 0; p < globalMesh.Cell0DTotalNumber(); p++)
      {
        dofsType[p] = (double)dofs.Cell0Ds_DOF[p].Type;
        dofsIndex[p] = dofs.Cell0Ds_DOF[p].GlobalIndex;
      }

      exporter.AddPoints(globalMesh.Cell0DsCoordinates(),
                         {
                           {
                             "DofType",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(dofsType.size()),
                             dofsType.data()
                           },
                           {
                             "DofIndex",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(dofsIndex.size()),
                             dofsIndex.data()
                           }
                         });

      exporter.Export(exportFolder +
                      "/DOFs_" +
                      "Mesh_Cell0Ds" +
                      ".vtu");
    }
  }
  // ***************************************************************************
  void DD_Utilities::Assemble(const int& rank,
                              const Problem_Info& problem_info,
                              const Gedim::IMeshDAO& globalMesh,
                              const std::vector<Eigen::MatrixXd>& squaresVertices,
                              const std::vector<double>& squaresArea,
                              const DOF_Info& dofs)
  {
    Eigen::MatrixXd referenceSquare_quadraturePoints;
    Eigen::VectorXd referenceSquare_quadratureWeights;
    Gedim::Quadrature_Gauss2D_Square::FillPointsAndWeights(2,
                                                           referenceSquare_quadraturePoints,
                                                           referenceSquare_quadratureWeights);

    Fem2DSquareLagrangePCC fem2D;
    const Fem2DSquareLagrangePCC::LocalSpace localSpace = fem2D.Compute();
    Eigen::MatrixXd referenceSquare_BasisFunctionValues = fem2D.Reference_BasisFunctions(localSpace,
                                                                                         referenceSquare_quadraturePoints);
    vector<Eigen::MatrixXd> referenceSquare_BasisFunctionDerivateValues = fem2D.Reference_BasisFunctionDerivatives(localSpace,
                                                                                                                   referenceSquare_quadraturePoints);

    SquareMapping mapping;

    for (unsigned int j_domain = 0; j_domain < problem_info.Num_1D_squares_domain; j_domain++)
    {
      for (unsigned int i_domain = 0; i_domain < problem_info.Num_1D_squares_domain; i_domain++)
      {
        // get mesh square
        const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                              i_domain,
                                                              j_domain,
                                                              problem_info.Num_1D_domains,
                                                              problem_info.Num_1D_squares_domain);

        const Square_Info square_info = Square_Info_Global(domain_info.Axes_Index.at(0),
                                                           domain_info.Axes_Index.at(1),
                                                           problem_info.Num_1D_points,
                                                           problem_info.Num_1D_squares);

        // get square geometric properties
        const Eigen::MatrixXd& square_Vertices = squaresVertices.at(square_info.Index);
        const double& square_Area = squaresArea.at(square_info.Index);

        // compute reference_element -> square map
        SquareMapping::Map map = mapping.Compute(square_Vertices,
                                                 square_Area);

        // map quadrature points and weights
        const Eigen::MatrixXd square_quadraturePoints = mapping.F(map,
                                                                  referenceSquare_quadraturePoints);
        const Eigen::VectorXd square_quadratureWeights = referenceSquare_quadratureWeights * abs(map.detQ);

        // map reference_element FEM basis functions
        const Eigen::MatrixXd square_basisFunctions_Values =
            fem2D.Map_BasisFunctions(localSpace,
                                     map,
                                     referenceSquare_BasisFunctionValues);

        const std::vector<Eigen::MatrixXd> square_basisFunctions_derivativeValues =
            fem2D.Map_BasisFunctionDerivatives(localSpace,
                                               map,
                                               referenceSquare_BasisFunctionDerivateValues);

        // compute equation elements
        const Eigen::VectorXd forcingTermValues = PDE_Equation::ForcingTerm(square_quadraturePoints);

        const Eigen::MatrixXd A = PDE_Equation::ComputeStiffnessMatrix(4,
                                                                       square_basisFunctions_derivativeValues,
                                                                       square_quadratureWeights);
        const Eigen::VectorXd f = PDE_Equation::ComputeCellForcingTerm(forcingTermValues,
                                                                       square_basisFunctions_Values,
                                                                       square_quadratureWeights);

      }
    }
  }
  // ***************************************************************************
}
