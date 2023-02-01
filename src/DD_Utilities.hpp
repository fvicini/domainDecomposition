#ifndef __DD_Utilities_H
#define __DD_Utilities_H

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
          unsigned int Num_2D_points;
          unsigned int Num_2D_squares;
          unsigned int Num_1D_points;
          unsigned int Num_1D_squares;
          unsigned int Num_1D_domains;
          unsigned int Num_1D_points_domain;
          unsigned int Num_1D_squares_domain;
      };

      struct Point_Info final
      {
          unsigned int Index;
      };

      struct Square_Info final
      {
          unsigned int Index;
          std::vector<unsigned int> Points_Index;
      };

      struct Domain_Info final
      {
          std::vector<unsigned int> Axes_Index;
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
              unsigned int GlobalIndex;
          };

          unsigned int Num_Dirichlets;
          unsigned int Num_Internals;
          unsigned int Num_Gamma;
          unsigned int Num_Globals;
          std::vector<DOF> Cell0Ds_DOF;
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

      static void ComputeErrors(const int& rank,
                                const Problem_Info& problem_info,
                                const Gedim::IMeshDAO& globalMesh,
                                const std::vector<Eigen::MatrixXd>& squaresVertices,
                                const std::vector<double>& squaresArea,
                                const DOF_Info& dofs,
                                const Gedim::IArray& internalSolution,
                                Eigen::VectorXd& errorL2,
                                Eigen::VectorXd& errorH1);

      static void ExportSolutionToVtu(const int& rank,
                                      const Problem_Info& problem_info,
                                      const DOF_Info& dofs,
                                      const Gedim::IMeshDAO& globalMesh,
                                      const Gedim::IArray& internalSolution,
                                      const std::string& exportFolder);

      static void ExportErrorToStream(const int& rank,
                                      const unsigned int& femOrder,
                                      const unsigned int& numCell2Ds,
                                      const unsigned int& numInternals,
                                      const unsigned int& numGamma,
                                      const double& h,
                                      const double& errorL2,
                                      const double& errorH1,
                                      const bool& printHeader,
                                      std::ostream& out,
                                      const char& separator);
  };

}

#endif
