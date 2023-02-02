#ifndef __DD_Utilities_H
#define __DD_Utilities_H

#include "ILinearSolver.hpp"
#include "ISparseArray.hpp"
#include "IArray.hpp"
#include "IMeshDAO.hpp"
#include "MeshUtilities.hpp"

#include <vector>

namespace DOMAIN_DECOMPOSITION
{
  class DD_Utilities final
  {
    public:
      struct Problem_Info
      {
          unsigned int Num_2D_points = 0;
          unsigned int Num_2D_squares = 0;
          unsigned int Num_1D_points = 0;
          unsigned int Num_1D_squares = 0;
          unsigned int Num_1D_domains = 0;
          unsigned int Num_1D_points_domain = 0;
          unsigned int Num_1D_squares_domain = 0;
      };

      struct Point_Info final
      {
          unsigned int Index = 0;
      };

      struct Square_Info final
      {
          unsigned int Index = 0;
          std::vector<unsigned int> Points_Index = {};
      };

      struct Domain_Info final
      {
          std::vector<unsigned int> Axes_Index = {};
      };

      struct DOF_Info final
      {
          struct DOF
          {
              enum struct Types
              {
                Unknwon = 0,
                Dirichlet = 1,
                Internal = 2,
                Gamma = 3
              };

              Types Type = Types::Unknwon;
              unsigned int Global_Index = 0; // index global in all domains
              unsigned int Local_Index = 0; // index local in the domain
          };

          struct Domain_DOF
          {
              unsigned int Num_Internals = 0;
              unsigned int Starting_Index = 0;
          };

          unsigned int Num_Dirichlets = 0;
          unsigned int Num_Internals = 0;
          unsigned int Num_Gamma = 0;
          unsigned int Num_Globals = 0;
          std::vector<Domain_DOF> Domains_DOF = {};
          std::vector<DOF> Cell0Ds_DOF = {};
      };

    private:
      /// \param i the x-axis index, from 0 to n_1D_points - 1
      /// \param j the y-axis index, from 0 to n_1D_points - 1
      /// \param n_1D_points the number of points on 1D
      /// \return the point information on global mesh
      static Point_Info Point_Info_Global(const unsigned int& i,
                                          const unsigned int& j,
                                          const unsigned int& n_1D_points);

      /// \param i the x-axis index, from 0 to n_1D_squares - 1
      /// \param j the y-axis index, from 0 to n_1D_squares - 1
      /// \param n_1D_points the number of points on 1D
      /// \param n_1D_squares the number of squares on 1D
      /// \return the square information on global mesh
      static Square_Info Square_Info_Global(const unsigned int& i,
                                            const unsigned int& j,
                                            const unsigned int& n_1D_points,
                                            const unsigned int& n_1D_squares);

      /// \param i_domain the x-axis index, from 0 to n_1D_squares_domain - 1
      /// \param j_domain the y-axis index, from 0 to n_1D_squares_domain - 1
      /// \param domainIndex the domain index, from 0 to n_domains
      /// \param n_1D_domains the number of domain on 1D
      /// \param n_1D_squares_domain the number of squares on 1D per domain
      /// \return the domain information on global mesh
      static Domain_Info GetDomain_Info_Global(const int& rank,
                                               const unsigned int& i_domain,
                                               const unsigned int& j_domain,
                                               const unsigned int& n_1D_domains,
                                               const unsigned int& n_1D_squares_domain);

      /// \brief Given an array p it computes the Schur complement product S*p
      static void ApplyShurToArray(const int& rank,
                                   const DOF_Info& dofs,
                                   const Gedim::ILinearSolver& A_II_solver,
                                   const Gedim::ISparseArray& A_IG,
                                   const Gedim::ISparseArray& A_GI,
                                   const Gedim::ISparseArray& A_GG,
                                   const Gedim::IArray& p,
                                   Gedim::IArray& Sp);

      /// \brief Solve S u_G = g with iterative algorithm
      static void ShurSolver(const int& rank,
                             const DOF_Info& dofs,
                             const Gedim::ILinearSolver& A_II_solver,
                             const Gedim::ISparseArray& A_IG,
                             const Gedim::ISparseArray& A_GI,
                             const Gedim::ISparseArray& A_GG,
                             const Gedim::IArray& g,
                             const unsigned int& max_iterations,
                             const double& tolerance,
                             const bool& conjugate,
                             Gedim::IArray& u_G);

      static void AII_Solver(const int& rank,
                             const DOF_Info& dofs,
                             const Gedim::ISparseArray& A_II,
                             const Gedim::ISparseArray& A_IG,
                             const Gedim::ISparseArray& A_GI,
                             const Gedim::ISparseArray& A_GG,
                             const Gedim::IArray& f_I,
                             const unsigned int& max_iterations,
                             const double& tolerance,
                             const bool& conjugate,
                             Gedim::IArray& u_G);

      static void PrintArray(const int& rank,
                             const std::string& v_name,
                             const Gedim::IArray& v);

    public:
      DD_Utilities();
      ~DD_Utilities();

      static void PrintMessage(const int& rank,
                               std::ostream& output,
                               const std::string& message,
                               const bool& onlyMaster);

      static void Assert(const int& rank,
                         const bool& result,
                         const std::string& message = "");

      static double StartTime();
      static double StopTime(const int& rank,
                             const int& n_domains,
                             const double& startTime,
                             const std::string& label,
                             const std::string& exportFolder);
      static void ExportTimes(const int& rank,
                              const int& n_domains,
                              const std::vector<double>& elapsedTimes,
                              const std::string& exportFolder);

      static Problem_Info ComputeProblemInfo(const int& rank,
                                             const int& n_domains,
                                             const Gedim::IMeshDAO& globalMesh);

      static void ExportDomainToVtu(const int& rank,
                                    const Problem_Info& problem_info,
                                    const Gedim::IMeshDAO& globalMesh,
                                    const std::string& exportFolder);

      static DOF_Info CreateDOFs(const int& rank,
                                 const int& n_domains,
                                 const Problem_Info& problem_info,
                                 const Gedim::IMeshDAO& globalMesh);

      static void ExportDOFsToVtu(const DOF_Info& dofs,
                                  const Gedim::IMeshDAO& globalMesh,
                                  const std::string& exportFolder);

      static void Assemble(const int& rank,
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
                           Gedim::IArray& f_G);

      static void Solve(const int& rank,
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
                        Gedim::IArray& u_G);

      static void ComputeErrors(const int& rank,
                                const Problem_Info& problem_info,
                                const Gedim::IMeshDAO& globalMesh,
                                const std::vector<Eigen::MatrixXd>& squaresVertices,
                                const std::vector<double>& squaresArea,
                                const DOF_Info& dofs,
                                const Gedim::IArray& u_I,
                                const Gedim::IArray& u_G,
                                Eigen::VectorXd& errorL2,
                                Eigen::VectorXd& errorH1);

      static void ExportSolutionToVtu(const int& rank,
                                      const Problem_Info& problem_info,
                                      const DOF_Info& dofs,
                                      const Gedim::IMeshDAO& globalMesh,
                                      const Gedim::IArray& u_I,
                                      const Gedim::IArray& u_G,
                                      const std::string& exportFolder);

      static void ExportErrorToStream(const int& rank,
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
                                      const char& separator);
  };

}

#endif
