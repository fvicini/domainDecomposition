#include "DD_Utilities.hpp"

#include "VTKUtilities.hpp"

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
  DD_Utilities::Square_Info DD_Utilities::GetSquare_Info_Global(const unsigned int& i,
                                                                const unsigned int& j,
                                                                const unsigned int& n_1D_points,
                                                                const unsigned int& n_1D_squares)
  {
    Square_Info info;

    info.Index = n_1D_squares * j + i;
    info.Points_Index.resize(4);
    info.Points_Index[0] = n_1D_points * j + i;
    info.Points_Index[1] = n_1D_points * j + (i + 1);
    info.Points_Index[2] = n_1D_points * (j + 1) + (i + 1);
    info.Points_Index[3] = n_1D_points * (j + 1) + i;

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
  void DD_Utilities::ExportDomainToVtu(const int& rank,
                                       const unsigned int& n_1D_points,
                                       const unsigned int& n_1D_squares,
                                       const unsigned int& n_1D_domains,
                                       const unsigned int& n_1D_squares_domain,
                                       const Gedim::IMeshDAO& globalMesh,
                                       const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      for (unsigned int i_domain = 0; i_domain < n_1D_squares_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < n_1D_squares_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                n_1D_domains,
                                                                n_1D_squares_domain);

          const Square_Info square_info = GetSquare_Info_Global(domain_info.Axes_Index.at(0),
                                                                domain_info.Axes_Index.at(1),
                                                                n_1D_points,
                                                                n_1D_squares);

          for (unsigned int v = 0; v < 4; v++)
          {
            const vector<double> cell0DDomain(1, rank);
            const vector<double> cell0DGlobalIndex(1, square_info.Points_Index[v]);
            exporter.AddPoint(globalMesh.Cell0DCoordinates(square_info.Points_Index[v]),
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
      }

      exporter.Export(exportFolder +
                      "/Domain_" +
                      to_string(rank) +
                      "_Mesh_Cell0Ds" +
                      ".vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      for (unsigned int i_domain = 0; i_domain < n_1D_squares_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < n_1D_squares_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                n_1D_domains,
                                                                n_1D_squares_domain);

          const Square_Info square_info = GetSquare_Info_Global(domain_info.Axes_Index.at(0),
                                                                domain_info.Axes_Index.at(1),
                                                                n_1D_points,
                                                                n_1D_squares);

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
}
