#ifndef __DD_Utilities_H
#define __DD_Utilities_H

#include "IMeshDAO.hpp"

#include <vector>

namespace DOMAIN_DECOMPOSITION
{
  class DD_Utilities final
  {
    public:
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

    private:

    public:
      DD_Utilities();
      ~DD_Utilities();

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

      static void ExportDomainToVtu(const int& rank,
                                    const unsigned int& n_1D_points,
                                    const unsigned int& n_1D_squares,
                                    const unsigned int& n_1D_domains,
                                    const unsigned int& n_1D_points_domain,
                                    const unsigned int& n_1D_squares_domain,
                                    const Gedim::IMeshDAO& globalMesh,
                                    const std::string& exportFolder);

      static void CreateDOFs(const unsigned int& n_1D_points,
                             const unsigned int& n_1D_squares,
                             const unsigned int& n_1D_domains,
                             const Gedim::IMeshDAO& globalMesh);
  };

}

#endif
