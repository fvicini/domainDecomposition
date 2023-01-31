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
  void DD_Utilities::ExportDomainToVtu(const int& rank,
                                       const unsigned int& n_1D_points,
                                       const unsigned int& n_1D_squares,
                                       const unsigned int& n_1D_domains,
                                       const unsigned int& n_1D_points_domain,
                                       const unsigned int& n_1D_squares_domain,
                                       const Gedim::IMeshDAO& globalMesh,
                                       const std::string& exportFolder)
  {
    {
      Gedim::VTKUtilities exporter;

      for (unsigned int i_domain = 0; i_domain < n_1D_points_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < n_1D_points_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                n_1D_domains,
                                                                n_1D_squares_domain);

          const Point_Info point_info = Point_Info_Global(domain_info.Axes_Index.at(0),
                                                          domain_info.Axes_Index.at(1),
                                                          n_1D_points);

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

      for (unsigned int i_domain = 0; i_domain < n_1D_squares_domain; i_domain++)
      {
        for (unsigned int j_domain = 0; j_domain < n_1D_squares_domain; j_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                n_1D_domains,
                                                                n_1D_squares_domain);

          const Square_Info square_info = Square_Info_Global(domain_info.Axes_Index.at(0),
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
  DD_Utilities::DOF_Info DD_Utilities::CreateDOFs(const int& n_domains,
                                                  const unsigned int& n_1D_points,
                                                  const unsigned int& n_1D_squares,
                                                  const unsigned int& n_1D_domains,
                                                  const unsigned int& n_1D_points_domain,
                                                  const unsigned int& n_1D_squares_domain,
                                                  const Gedim::IMeshDAO& globalMesh)
  {
    DOF_Info info;

    // initialize
    info.Num_Dirichlets = 0; // all external borders are dirichlets
    info.Num_Internals = 0;
    info.Num_Gamma = 0;
    info.Num_Globals = globalMesh.Cell0DTotalNumber();
    info.Cell0Ds_Type.resize(globalMesh.Cell0DTotalNumber(), DOF_Info::Types::Unknwon);
    info.Cell0Ds_GlobalIndex.resize(globalMesh.Cell0DTotalNumber(), 0);

    // check mesh points type
    for (unsigned int i = 0; i < n_1D_points; i++)
    {
      for (unsigned int j = 0; j < n_1D_points; j++)
      {
        const Point_Info global_info = Point_Info_Global(i, j, n_1D_points);

        if (i == 0 || i == n_1D_points - 1 ||
            j == 0 || j == n_1D_points - 1)
        {
          info.Cell0Ds_Type[global_info.Index] = DOF_Info::Types::Dirichlet;
          info.Num_Dirichlets++;
        }
        else if (i % n_1D_squares_domain == 0 ||
                 j % n_1D_squares_domain == 0)
        {
          info.Cell0Ds_Type[global_info.Index] = DOF_Info::Types::Gamma;
          info.Num_Gamma++;
        }
        else
        {
          info.Cell0Ds_Type[global_info.Index] = DOF_Info::Types::Internal;
          info.Num_Internals++;
        }
      }
    }

    // numerate dirichlet and gamma points
    unsigned int dirichlet_counter = 0;
    unsigned int gamma_counter = 0;
    for (unsigned int p = 0; p < globalMesh.Cell0DTotalNumber(); p++)
    {
      switch (info.Cell0Ds_Type[p])
      {
        case DOF_Info::Types::Dirichlet:
          info.Cell0Ds_GlobalIndex[p] = info.Num_Internals +
                                        info.Num_Gamma +
                                        dirichlet_counter++;
          break;
        case DOF_Info::Types::Gamma:
          info.Cell0Ds_GlobalIndex[p] = info.Num_Internals +
                                        gamma_counter++;
          break;
        case DOF_Info::Types::Internal:
          continue;
        default:
          throw runtime_error("unkwnon point type");
      }
    }

    Gedim::Output::Assert(dirichlet_counter == info.Num_Dirichlets);
    Gedim::Output::Assert(gamma_counter == info.Num_Gamma);

    // numerate internal dofs with processes order
    unsigned int internal_counter = 0;
    for (unsigned int rank = 0; rank < n_domains; rank++)
    {
      for (unsigned int j_domain = 0; j_domain < n_1D_points_domain; j_domain++)
      {
        for (unsigned int i_domain = 0; i_domain < n_1D_points_domain; i_domain++)
        {
          const Domain_Info domain_info = GetDomain_Info_Global(rank,
                                                                i_domain,
                                                                j_domain,
                                                                n_1D_domains,
                                                                n_1D_squares_domain);

          const Point_Info point_info = Point_Info_Global(domain_info.Axes_Index.at(0),
                                                          domain_info.Axes_Index.at(1),
                                                          n_1D_points);

          switch (info.Cell0Ds_Type[point_info.Index])
          {
            case DOF_Info::Types::Dirichlet:
            case DOF_Info::Types::Gamma:
              continue;
            case DOF_Info::Types::Internal:
              info.Cell0Ds_GlobalIndex[point_info.Index] = internal_counter++;
              break;
            default:
              throw runtime_error("unkwnon point type");
          }
        }
      }
    }

    Gedim::Output::Assert(internal_counter == info.Num_Internals);

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
        dofsType[p] = (double)dofs.Cell0Ds_Type[p];
        dofsIndex[p] = dofs.Cell0Ds_GlobalIndex[p];
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
}
