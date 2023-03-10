#include "DD_Utilities.hpp"

#include <numeric>

#include "Eigen_CholeskySolver.hpp"
#include "Eigen_LUSolver.hpp"
#include "Eigen_Array.hpp"
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
  double DD_Utilities::StartTime()
  {
    return MPI::Wtime();
  }
  // ***************************************************************************
  double DD_Utilities::StopTime(const int& rank,
                                const int& n_domains,
                                const double& startTime,
                                const std::string& label,
                                const std::string& exportFolder)
  {
    double elapsedTime = MPI::Wtime() - startTime;
    double maxElapsedTime = 0.0;
    MPI_Reduce(&elapsedTime, &maxElapsedTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
      // export max elapsed time
      string filePath = exportFolder + "/Time.csv";
      const bool fileExists = Gedim::Output::FileExists(filePath);
      ofstream outFile(filePath,
                       ios_base::app);
      outFile.precision(16);

      if(!fileExists)
      {
        outFile<< "Numdomains"<< " ";
        outFile<< "Label"<< " ";
        outFile<< "ElapsedTime"<< endl;
      }

      outFile<< scientific<< n_domains<< " ";
      outFile<< scientific<< label<< " ";
      outFile<< scientific<< maxElapsedTime<< endl;
    }

    {
      const string processesTimeFolder = exportFolder + "/Processes";
      Gedim::Output::CreateFolder(processesTimeFolder);
      // export local elapsed time
      string filePath = processesTimeFolder + "/Time_" + to_string(rank) + ".csv";
      const bool fileExists = Gedim::Output::FileExists(filePath);
      ofstream outFile(filePath,
                       ios_base::app);
      outFile.precision(16);

      if(!fileExists)
      {
        outFile<< "Numdomains"<< " ";
        outFile<< "Label"<< " ";
        outFile<< "ElapsedTime"<< endl;
      }

      outFile<< scientific<< n_domains<< " ";
      outFile<< scientific<< label<< " ";
      outFile<< scientific<< elapsedTime<< endl;
    }

    return elapsedTime;
  }
  // ***************************************************************************
  void DD_Utilities::ExportTimes(const int& rank,
                                 const int& n_domains,
                                 const std::vector<double>& elapsedTimes,
                                 const std::string& exportFolder)
  {

    const double totalTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0);
    double maxTotalTime = 0.0;
    MPI_Reduce(&totalTime, &maxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
      // export max total elapsed time
      string filePath = exportFolder + "/Time.csv";
      const bool fileExists = Gedim::Output::FileExists(filePath);
      ofstream outFile(filePath,
                       ios_base::app);
      outFile.precision(16);

      if(!fileExists)
      {
        outFile<< "Numdomains"<< " ";
        outFile<< "Label"<< " ";
        outFile<< "ElapsedTime"<< endl;
      }

      outFile<< scientific<< n_domains<< " ";
      outFile<< scientific<< "Total"<< " ";
      outFile<< scientific<< maxTotalTime<< endl;
    }

    {
      // export local elapsed time
      const string processesTimeFolder = exportFolder + "/Processes";
      Gedim::Output::CreateFolder(processesTimeFolder);
      string filePath = processesTimeFolder + "/Time_" + to_string(rank) + ".csv";
      const bool fileExists = Gedim::Output::FileExists(filePath);
      ofstream outFile(filePath,
                       ios_base::app);
      outFile.precision(16);

      if(!fileExists)
      {
        outFile<< "Numdomains"<< " ";
        outFile<< "Label"<< " ";
        outFile<< "ElapsedTime"<< endl;
      }

      outFile<< scientific<< n_domains<< " ";
      outFile<< scientific<< "Total"<< " ";
      outFile<< scientific<< totalTime<< endl;
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
  void DD_Utilities::ApplyShurToArray(const int& rank,
                                      const DOF_Info& dofs,
                                      const Gedim::ILinearSolver& A_II_solver,
                                      const Gedim::ISparseArray& A_IG,
                                      const Gedim::ISparseArray& A_GI,
                                      const Gedim::ISparseArray& A_GG,
                                      const Gedim::IArray& p,
                                      Gedim::IArray& Sp)
  {
    Gedim::Eigen_Array<> w_I, rhs_I;

    w_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);
    rhs_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);
    Sp.SetSize(dofs.Num_Gamma);

    // solve A_II w_I = - A_IG * p
    rhs_I.SubtractionMultiplication(A_IG, p);
    A_II_solver.Solve(rhs_I, w_I);

    // compute Sp = sum_domain A_GI * w_I + A_GG * p
    Sp.SumMultiplication(A_GI, w_I);
    Sp.SumMultiplication(A_GG, p);
    MPI_Allreduce(MPI_IN_PLACE, Sp.Data(), Sp.Size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  // ***************************************************************************
  void DD_Utilities::ShurSolver(const int& rank,
                                const DOF_Info& dofs,
                                const Gedim::ILinearSolver& A_II_solver,
                                const Gedim::ISparseArray& A_IG,
                                const Gedim::ISparseArray& A_GI,
                                const Gedim::ISparseArray& A_GG,
                                const Gedim::IArray& g,
                                const unsigned int& max_iterations,
                                const double& tolerance,
                                const bool& conjugate,
                                Gedim::IArray& u_G)
  {
    Gedim::Eigen_Array<> r_k, p_k, Sp_k;
    r_k.Copy(g);
    p_k.Copy(r_k);

    double beta_k = 0.0, alpha_k = 0.0;
    unsigned int iteration = 0;

    double r_0_dot = r_k.Dot(r_k);
    double r_k_1_dot = r_0_dot;
    double r_k_dot = r_0_dot;

    if (rank == 0)
    {
      cout<< scientific<< "Schur Solver"<< " ";
      cout<< "it "<< iteration<< "/"<< max_iterations<< " ";
      cout<< "res "<< sqrt(r_k_dot)<< "/"<< tolerance * sqrt(r_0_dot)<< endl;
    }

    while (r_k_dot > tolerance * tolerance * r_0_dot &&
           iteration < max_iterations)
    {
      beta_k = (!conjugate || iteration == 0) ? 0.0 : r_k_dot / r_k_1_dot;

      // compute p_k = r_k + beta_k * p_k_1
      p_k *= beta_k;
      p_k += r_k; // Conjugate Gradient

      // compute S * p_k
      ApplyShurToArray(rank,
                       dofs,
                       A_II_solver,
                       A_IG,
                       A_GI,
                       A_GG,
                       p_k,
                       Sp_k);

      alpha_k = r_k_dot / p_k.Dot(Sp_k);

      // compute u_G_k = u_G_k_1 + alpha_k * p_k
      u_G += (p_k * alpha_k);

      // compute r_k = r_k_1 - alpha_k * S * p_k
      r_k_1_dot = r_k_dot;

      r_k -= (Sp_k * alpha_k);
      r_k_dot = r_k.Dot(r_k);

      iteration++;

      if (rank == 0)
      {
        cout<< scientific<< "Schur Solver"<< " ";
        cout<< "it "<< iteration<< "/"<< max_iterations<< " ";
        cout<< "res "<< sqrt(r_k_dot)<< "/"<< tolerance * sqrt(r_0_dot)<< endl;
      }
    }
  }
  // ***************************************************************************
  void DD_Utilities::AII_Solver(const int& rank,
                                const DOF_Info& dofs,
                                const Gedim::ISparseArray& A_II,
                                const Gedim::ISparseArray& A_IG,
                                const Gedim::ISparseArray& A_GI,
                                const Gedim::ISparseArray& A_GG,
                                const Gedim::IArray& f_I,
                                const unsigned int& max_iterations,
                                const double& tolerance,
                                const bool& conjugate,
                                Gedim::IArray& u_I)
  {
    Gedim::Eigen_Array<> r_k, p_k;
    r_k.Copy(f_I);
    p_k.Copy(r_k);

    Gedim::Eigen_Array<> A_pk;
    A_pk.SetSize(dofs.Domains_DOF[rank].Num_Internals);

    double beta_k = 0.0, alpha_k = 0.0;
    unsigned int iteration = 0;

    double r_0_dot = r_k.Dot(r_k);
    double r_k_1_dot = r_0_dot;
    double r_k_dot = r_0_dot;

    while (r_k_dot > tolerance * tolerance * r_0_dot &&
           iteration < max_iterations)
    {
      if (rank == 0)
      {
        cout<< scientific<< "AII Solver"<< " ";
        cout<< "it "<< iteration<< "/"<< max_iterations<< " ";
        cout<< "res "<< sqrt(r_k_dot)<< "/"<< tolerance * sqrt(r_0_dot)<< endl;
      }

      // compute beta_k = <r_k, r_k> / <r_k_1, r_k_1>
      beta_k = (!conjugate || iteration == 0) ? 0.0 : r_k_dot / r_k_1_dot;

      // compute p_k = r_k + beta_k * p_k_1
      p_k *= beta_k;
      p_k += r_k; // Conjugate Gradient

      // compute alpha_k = <r_k, r_k> / <p_k, A * p_k>
      A_pk.Zeros();
      A_pk.SumMultiplication(A_II, p_k);
      alpha_k = r_k_dot / p_k.Dot(A_pk);

      // compute u_G_k = u_G_k_1 + alpha_k * p_k
      u_I += (p_k * alpha_k);

      // save <r_k_1, r_k_1>
      r_k_1_dot = r_k_dot;

      // compute r_k = r_k_1 - alpha_k * S * p_k
      r_k -= (A_pk * alpha_k);
      r_k_dot = r_k.Dot(r_k);

      iteration++;
    }
  }
  // ***************************************************************************
  void DD_Utilities::PrintArray(const int& rank,
                                const std::string& v_name,
                                const Gedim::IArray& v)
  {
    {
      using namespace Gedim;
      cout<< scientific<< "Process "<< rank<< " " + v_name + " "<< v<< endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
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
    A_II.SetSize(dofs.Domains_DOF[rank].Num_Internals,
                 dofs.Domains_DOF[rank].Num_Internals);
    A_IG.SetSize(dofs.Domains_DOF[rank].Num_Internals,
                 dofs.Num_Gamma);
    A_GI.SetSize(dofs.Num_Gamma,
                 dofs.Domains_DOF[rank].Num_Internals);
    A_GG.SetSize(dofs.Num_Gamma,
                 dofs.Num_Gamma);

    f_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);
    f_G.SetSize(dofs.Num_Gamma);

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
                           const unsigned int& max_iterations,
                           const double& tolerance,
                           const bool& conjugate,
                           Gedim::IArray& u_I,
                           Gedim::IArray& u_G)
  {
    u_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);
    u_G.SetSize(dofs.Num_Gamma);

    // compute cholesky factorization of A_II
    Gedim::Eigen_LUSolver<> A_II_solver;
    A_II_solver.Initialize(A_II);

    if (dofs.Num_Gamma == 0)
    {
      // solve A_II * u_I = f_I
      A_II_solver.Solve(f_I, u_I);
      return;
    }

    // initialize auxiliary variables
    Gedim::Eigen_Array<> h_I, g, rhs_I;

    h_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);
    g.SetSize(dofs.Num_Gamma);
    rhs_I.SetSize(dofs.Domains_DOF[rank].Num_Internals);

    // solve internal system A_II * h_I = f_I
    A_II_solver.Solve(f_I, h_I);

    // compute g = sum_domain f_G - A_GI * h_I
    g += f_G;
    g.SubtractionMultiplication(A_GI, h_I);
    MPI_Allreduce(MPI_IN_PLACE, g.Data(), g.Size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // solve S u_G = g with CG
    ShurSolver(rank,
               dofs,
               A_II_solver,
               A_IG,
               A_GI,
               A_GG,
               g,
               max_iterations,
               tolerance,
               conjugate,
               u_G);

    // solve A_II * u_i = f_I - A_IG * u_g
    rhs_I += f_I;
    rhs_I.SubtractionMultiplication(A_IG, u_G);

    A_II_solver.Solve(rhs_I, u_I);
  }
  // ***************************************************************************
  void DD_Utilities::ComputeErrors(const int& rank,
                                   const Problem_Info& problem_info,
                                   const Gedim::IMeshDAO& globalMesh,
                                   const std::vector<Eigen::MatrixXd>& squaresVertices,
                                   const std::vector<double>& squaresArea,
                                   const DOF_Info& dofs,
                                   const Gedim::IArray& u_I,
                                   const Gedim::IArray& u_G,
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

          if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
              DOF_Info::DOF::Types::Internal)
          {
            const unsigned int i_dom_loc = dofs.Cell0Ds_DOF[i_mesh_point_index].Local_Index;
            localNumericSolution[i_loc] = u_I[i_dom_loc];
          }
          else if (dofs.Cell0Ds_DOF[i_mesh_point_index].Type ==
                   DOF_Info::DOF::Types::Gamma)
          {
            const unsigned int i_glb = dofs.Cell0Ds_DOF[i_mesh_point_index].Global_Index;
            localNumericSolution[i_loc] = u_G[i_glb];
          }
          else
            throw runtime_error("Unknown J DOF type");
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
                                         const Gedim::IArray& u_I,
                                         const Gedim::IArray& u_G,
                                         const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      const unsigned int numPoints = problem_info.Num_1D_points_domain *
                                     problem_info.Num_1D_points_domain;

      Eigen::MatrixXd coordinates = Eigen::MatrixXd(3,
                                                    numPoints);

      vector<double> numericalSolution(numPoints, 0.0);

      for (unsigned int j_domain = 0; j_domain < problem_info.Num_1D_points_domain; j_domain++)
      {
        for (unsigned int i_domain = 0; i_domain < problem_info.Num_1D_points_domain; i_domain++)
        {
          const unsigned int point_index = j_domain * problem_info.Num_1D_points_domain + i_domain;

          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                problem_info.Num_1D_domains,
                                                                problem_info.Num_1D_squares_domain);

          const Point_Info point_info = Point_Info_Global(domain_info.Axes_Index.at(0),
                                                          domain_info.Axes_Index.at(1),
                                                          problem_info.Num_1D_points);

          coordinates.col(point_index)<< globalMesh.Cell0DCoordinates(point_info.Index);

          switch (dofs.Cell0Ds_DOF[point_info.Index].Type)
          {
            case DOF_Info::DOF::Types::Dirichlet:
              continue;
            case DOF_Info::DOF::Types::Gamma:
              numericalSolution[point_index] = u_G[
                                               dofs.Cell0Ds_DOF[point_info.Index].Global_Index];
              break;
            case DOF_Info::DOF::Types::Internal:
              numericalSolution[point_index] = u_I[
                                               dofs.Cell0Ds_DOF[point_info.Index].Local_Index];
              break;
            default:
              throw runtime_error("unkwnon point type");
          }
        }
      }

      const Eigen::VectorXd exactSolution = PDE_Equation::ExactSolution(coordinates);

      exporter.AddPoints(coordinates,
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
                      "_Mesh_Cell0Ds" +
                      ".vtu");
    }
  }
  // ***************************************************************************
  void DD_Utilities::ExportErrorToStream(const int& rank,
                                         const int& n_domains,
                                         const unsigned int& femOrder,
                                         const unsigned int& numCell2Ds,
                                         const unsigned int& numInternals,
                                         const unsigned int& numGamma,
                                         const double& h,
                                         const double& H,
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
      out<< "Numdomains" << separator;
      out<< "FemOrder" << separator;
      out<< "Cell2Ds" <<  separator;
      out<< "NumInternals" <<  separator;
      out<< "NumGamma" <<  separator;
      out<< "h" <<  separator;
      out<< "H" <<  separator;
      out<< "L2" <<  separator;
      out<< "H1" << endl;
    }

    out.precision(16);
    out<< scientific<< n_domains<< separator;
    out<< scientific<< femOrder<< separator;
    out<< scientific<< numCell2Ds<< separator;
    out<< scientific<< numInternals<< separator;
    out<< scientific<< numGamma<< separator;
    out<< scientific<< h << separator;
    out<< scientific<< H << separator;
    out<< scientific<< errorL2<< separator;
    out<< scientific<< errorH1<< endl;
  }
  // ***************************************************************************
}
