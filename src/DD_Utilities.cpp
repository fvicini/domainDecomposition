#include "DD_Utilities.hpp"

#include "Eigen_CholeskySolver.hpp"
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
    info.Domains_DOF.resize(n_domains);
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
          info.Cell0Ds_DOF[p].Global_Index = dirichlet_counter++;
          info.Cell0Ds_DOF[p].Local_Index = info.Cell0Ds_DOF[p].Global_Index;
          break;
        case DOF_Info::DOF::Types::Gamma:
          info.Cell0Ds_DOF[p].Global_Index = gamma_counter++;
          info.Cell0Ds_DOF[p].Local_Index = info.Cell0Ds_DOF[p].Global_Index;
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
      info.Domains_DOF[rank].Num_Internals = 0;
      info.Domains_DOF[rank].Starting_Index = internal_counter;
      unsigned int local_internal_counter = 0;

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
              info.Domains_DOF[rank].Num_Internals++;
              info.Cell0Ds_DOF[point_info.Index].Global_Index = internal_counter++;
              info.Cell0Ds_DOF[point_info.Index].Local_Index = local_internal_counter++;
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
      vector<double> dofsGlobalIndex(globalMesh.Cell0DTotalNumber());
      vector<double> dofsLocalIndex(globalMesh.Cell0DTotalNumber());

      for (unsigned int p = 0; p < globalMesh.Cell0DTotalNumber(); p++)
      {
        dofsType[p] = (double)dofs.Cell0Ds_DOF[p].Type;
        dofsGlobalIndex[p] = dofs.Cell0Ds_DOF[p].Global_Index;
        dofsLocalIndex[p] = dofs.Cell0Ds_DOF[p].Local_Index;
      }

      exporter.AddPoints(globalMesh.Cell0DsCoordinates(),
                         {
                           {
                             "Type",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(dofsType.size()),
                             dofsType.data()
                           },
                           {
                             "Global_Index",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(dofsGlobalIndex.size()),
                             dofsGlobalIndex.data()
                           },
                           {
                             "Local_Index",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(dofsLocalIndex.size()),
                             dofsLocalIndex.data()
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
                              const DOF_Info& dofs,
                              Gedim::ISparseArray& A_II,
                              Gedim::ISparseArray& A_IG,
                              Gedim::ISparseArray& A_GI,
                              Gedim::ISparseArray& A_GG,
                              Gedim::IArray& f_I,
                              Gedim::IArray& f_G)
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

        const Eigen::MatrixXd A = PDE_Equation::ComputeStiffnessMatrix(localSpace.NumberBasisFunctions,
                                                                       square_basisFunctions_derivativeValues,
                                                                       square_quadratureWeights);
        const Eigen::VectorXd f = PDE_Equation::ComputeCellForcingTerm(forcingTermValues,
                                                                       square_basisFunctions_Values,
                                                                       square_quadratureWeights);

        // insert local matrix into global matrix
        for (unsigned int i_loc = 0; i_loc < localSpace.NumberBasisFunctions; i_loc++)
        {
          const unsigned int i_mesh_point_index = square_info.Points_Index.at(i_loc);

          if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
              DOF_Info::DOF::Types::Dirichlet)
            continue;

          if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
              DOF_Info::DOF::Types::Internal)
          {
            const unsigned int i_dom_loc = dofs.Cell0Ds_DOF[i_mesh_point_index].Local_Index;

            f_I.AddValue(i_dom_loc,
                         f[i_loc]);

            for (unsigned int j_loc = 0; j_loc < localSpace.NumberBasisFunctions; j_loc++)
            {
              const unsigned int j_mesh_point_index = square_info.Points_Index.at(j_loc);

              if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                  DOF_Info::DOF::Types::Dirichlet)
                continue;

              if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                  DOF_Info::DOF::Types::Internal)
              {
                const unsigned int j_dom_loc = dofs.Cell0Ds_DOF[j_mesh_point_index].Local_Index;

                A_II.Triplet(i_dom_loc,
                             j_dom_loc,
                             A(i_loc, j_loc));
              }
              else if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                       DOF_Info::DOF::Types::Gamma)
              {
                const unsigned int j_glb = dofs.Cell0Ds_DOF[j_mesh_point_index].Global_Index;

                A_IG.Triplet(i_dom_loc,
                             j_glb,
                             A(i_loc, j_loc));
              }
              else
                throw runtime_error("Unknown J DOF type");
            }
          }
          else if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
                   DOF_Info::DOF::Types::Gamma)
          {
            const unsigned int i_glb = dofs.Cell0Ds_DOF[i_mesh_point_index].Global_Index;

            f_G.AddValue(i_glb,
                         f[i_loc]);

            for (unsigned int j_loc = 0; j_loc < localSpace.NumberBasisFunctions; j_loc++)
            {
              const unsigned int j_mesh_point_index = square_info.Points_Index.at(j_loc);

              if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                  DOF_Info::DOF::Types::Dirichlet)
                continue;

              if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                  DOF_Info::DOF::Types::Internal)
              {
                const unsigned int j_dom_loc = dofs.Cell0Ds_DOF[j_mesh_point_index].Local_Index;

                A_GI.Triplet(i_glb,
                             j_dom_loc,
                             A(i_loc, j_loc));
              }
              else if (dofs.Cell0Ds_DOF[j_mesh_point_index].Type ==
                       DOF_Info::DOF::Types::Gamma)
              {
                const unsigned int j_glb = dofs.Cell0Ds_DOF[j_mesh_point_index].Global_Index;

                A_GG.Triplet(i_glb,
                             j_glb,
                             A(i_loc, j_loc));
              }
              else
                throw runtime_error("Unknown J DOF type");
            }
          }
          else
            throw runtime_error("Unknown I DOF type");
        }
      }
    }

    f_I.Create();
    f_G.Create();
    A_II.Create();
    A_IG.Create();
    A_GI.Create();
    A_GG.Create();
  }
  // ***************************************************************************
  void DD_Utilities::Solve(const int& rank,
                           const Problem_Info& problem_info,
                           const Gedim::IMeshDAO& globalMesh,
                           const DOF_Info& dofs,
                           const Gedim::ISparseArray& A_II,
                           const Gedim::ISparseArray& A_IG,
                           const Gedim::ISparseArray& A_GI,
                           const Gedim::ISparseArray& A_GG,
                           const Gedim::IArray& f_I,
                           const Gedim::IArray& f_G,
                           Gedim::IArray& u_I,
                           Gedim::IArray& u_G)
  {
    // solve internal system
    Gedim::Eigen_CholeskySolver<> choleskySolver;
    choleskySolver.Initialize(A_II,
                              f_I,
                              u_I);
    choleskySolver.Solve();
  }
  // ***************************************************************************
  void DD_Utilities::ComputeErrors(const int& rank,
                                   const Problem_Info& problem_info,
                                   const Gedim::IMeshDAO& globalMesh,
                                   const std::vector<Eigen::MatrixXd>& squaresVertices,
                                   const std::vector<double>& squaresArea,
                                   const DOF_Info& dofs,
                                   const Gedim::IArray& internalSolution,
                                   Eigen::VectorXd& errorL2,
                                   Eigen::VectorXd& errorH1)
  {
    errorL2 = Eigen::VectorXd::Zero(globalMesh.Cell2DTotalNumber());
    errorH1 = Eigen::VectorXd::Zero(globalMesh.Cell2DTotalNumber());

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
        const Eigen::VectorXd exactSolution = PDE_Equation::ExactSolution(square_quadraturePoints);
        const Eigen::VectorXd exactDerivativeSolution_X = PDE_Equation::ExactDerivativeSolution(0,
                                                                                                square_quadraturePoints);
        const Eigen::VectorXd exactDerivativeSolution_Y = PDE_Equation::ExactDerivativeSolution(1,
                                                                                                square_quadraturePoints);

        // compute local numerical solution
        Eigen::VectorXd localNumericSolution = Eigen::VectorXd::Zero(localSpace.NumberBasisFunctions);
        for (unsigned int i_loc = 0; i_loc < localSpace.NumberBasisFunctions; i_loc++)
        {
          const unsigned int i_mesh_point_index = square_info.Points_Index.at(i_loc);

          if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
              DOF_Info::DOF::Types::Dirichlet)
            continue;

          const unsigned int i_glb = dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
                                     DOF_Info::DOF::Types::Internal ?
                                       dofs.Cell0Ds_DOF[i_mesh_point_index].Global_Index :
                                       dofs.Num_Internals +
                                       dofs.Cell0Ds_DOF[i_mesh_point_index].Global_Index;

          localNumericSolution[i_loc] = internalSolution[i_glb];
        }

        Eigen::VectorXd localErrorL2 = (square_basisFunctions_Values * localNumericSolution -
                                        exactSolution).array().square();

        Eigen::VectorXd localErrorH1 = Eigen::VectorXd::Zero(square_quadraturePoints.cols());
        localErrorH1.array() += (square_basisFunctions_derivativeValues[0] * localNumericSolution -
                                exactDerivativeSolution_X).array().square();
        localErrorH1.array() += (square_basisFunctions_derivativeValues[1] * localNumericSolution -
                                exactDerivativeSolution_Y).array().square();

        errorL2[square_info.Index] = square_quadratureWeights.transpose() * localErrorL2;
        errorH1[square_info.Index] = square_quadratureWeights.transpose() * localErrorH1;
      }
    }
  }
  // ***************************************************************************
  void DD_Utilities::ExportSolutionToVtu(const int& rank,
                                         const Problem_Info& problem_info,
                                         const DOF_Info& dofs,
                                         const Gedim::IMeshDAO& globalMesh,
                                         const Gedim::IArray& internalSolution,
                                         const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      const Eigen::MatrixXd coordinates = globalMesh.Cell0DsCoordinates();
      Eigen::VectorXd exactSolution = PDE_Equation::ExactSolution(coordinates);
      vector<double> numericalSolution(globalMesh.Cell0DTotalNumber(), 0.0);

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

          switch (dofs.Cell0Ds_DOF[point_info.Index].Type)
          {
            case DOF_Info::DOF::Types::Dirichlet:
              continue;
            case DOF_Info::DOF::Types::Gamma:
              numericalSolution[point_info.Index] = internalSolution[
                                                    dofs.Num_Internals +
                                                    dofs.Cell0Ds_DOF[point_info.Index].Global_Index];
              break;
            case DOF_Info::DOF::Types::Internal:
              numericalSolution[point_info.Index] = internalSolution[
                                                    dofs.Cell0Ds_DOF[point_info.Index].Global_Index];
              break;
            default:
              throw runtime_error("unkwnon point type");
          }

        }
      }

      exporter.AddPolygons(coordinates,
                           globalMesh.Cell2DsVertices(),
                           {
                             {
                               "Numerical",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(numericalSolution.size()),
                               numericalSolution.data()
                             },
                             {
                               "Exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(exactSolution.size()),
                               exactSolution.data()
                             }
                           });

      exporter.Export(exportFolder +
                      "/Solution_" +
                      to_string(rank) +
                      "_Mesh_Cell2Ds" +
                      ".vtu");
    }
  }
  // ***************************************************************************
  void DD_Utilities::ExportErrorToStream(const int& rank,
                                         const unsigned int& femOrder,
                                         const unsigned int& numCell2Ds,
                                         const unsigned int& numInternals,
                                         const unsigned int& numGamma,
                                         const double& h,
                                         const double& errorL2,
                                         const double& errorH1,
                                         const bool& printHeader,
                                         std::ostream& out,
                                         const char& separator)
  {
    if (rank > 0)
      return;

    if (printHeader)
    {
      out<< "FemOrder" << separator;
      out<< "Cell2Ds" <<  separator;
      out<< "NumInternals" <<  separator;
      out<< "NumGamma" <<  separator;
      out<< "h" <<  separator;
      out<< "L2" <<  separator;
      out<< "H1" << endl;
    }

    out.precision(16);
    out<< scientific<< femOrder<< separator;
    out<< scientific<< numCell2Ds<< separator;
    out<< scientific<< numInternals<< separator;
    out<< scientific<< numGamma<< separator;
    out<< scientific<< h << separator;
    out<< scientific<< errorL2<< separator;
    out<< scientific<< errorH1<< endl;
  }
  // ***************************************************************************
}
