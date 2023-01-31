#include "EllipticProblem.hpp"

#include "Configurations.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_CholeskySolver.hpp"

#include "mpi.h"

#include "DD_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace DOMAIN_DECOMPOSITION
{
  // ***************************************************************************
  EllipticProblem_ProgramConfiguration::EllipticProblem_ProgramConfiguration()
  {
    // Export parameters
    Gedim::Configurations::AddProperty("ExportFolder",
                                       "./Run",
                                       "Folder where to export data (Default: ./Export)");
    // Geometric parameters
    Gedim::Configurations::AddProperty("GeometricTolerance",
                                       1.0e-8,
                                       "Geometric tolerance to perform 1D operations (Default: machine epsilon)");

    // Mesh parameters
    Gedim::Configurations::AddProperty("MeshParameter",
                                       0.25,
                                       "Mesh Parameter (Default: 0.25)");
  }
  // ***************************************************************************
  void EllipticProblem::PrintMessage(const int& rank,
                                     ostream& output,
                                     const string& message,
                                     const bool& onlyMaster)
  {
    if (onlyMaster && rank != 0)
      return;

    output<< ">> "<< "Process "<< rank<< ": "<< message<< endl;
  }
  // ***************************************************************************
  EllipticProblem::EllipticProblem(const EllipticProblem_ProgramConfiguration& config) :
    config(config)
  {
  }
  EllipticProblem::~EllipticProblem()
  {
  }
  // ***************************************************************************
  void EllipticProblem::Run(const int& rank,
                            const int& n_domains) const
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = config.GeometricTolerance();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    /// Create folders˝
    const string exportFolder = config.ExportFolder();
    if (rank == 0)
      Gedim::Output::CreateFolder(exportFolder);

    const string exportCsvFolder = exportFolder + "/Mesh";
    if (rank == 0)
      Gedim::Output::CreateFolder(exportCsvFolder);
    const string exportVtuFolder = exportFolder + "/Paraview";
    if (rank == 0)
      Gedim::Output::CreateFolder(exportVtuFolder);
    const string exportVtuDomainFolder = exportVtuFolder + "/Domains";
    if (rank == 0)
      Gedim::Output::CreateFolder(exportVtuDomainFolder);
    const string exportSolutionFolder = exportFolder + "/Solution";
    if (rank == 0)
      Gedim::Output::CreateFolder(exportSolutionFolder);

    const string logFolder = exportFolder + "/Log";

    /// Export Configuration of the following Run
    if (rank == 0)
      Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini",
                                         false);

    /// Set Log folder
    if (rank == 0)
      Gedim::Output::CreateFolder(logFolder);
    Gedim::LogFile::LogFolder = logFolder;

    /// Create problem
    PrintMessage(rank, cout, "Create Domain...", false);

    const Eigen::Vector3d rectangleOrigin(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBase(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeight(0.0, 1.0, 0.0);

    const Eigen::MatrixXd domain = geometryUtilities.CreateParallelogram(rectangleOrigin,
                                                                         rectangleBase,
                                                                         rectangleHeight);

    PrintMessage(rank, cout, "Create Domain SUCCESS", false);

    // Export domain
    PrintMessage(rank, cout, "Export Domain...", true);
    if (rank == 0)
    {
      Gedim::VTKUtilities vtkUtilities;
      vtkUtilities.AddPolygon(domain);
      vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }
    PrintMessage(rank, cout, "Export Domain SUCCESS", true);

    /// Create domain mesh
    PrintMessage(rank, cout,
                 "Create Domain Mesh with parameter " +
                 to_string(config.MeshParameter()) +
                 "...", false);

    Gedim::MeshMatrices domainMeshData;
    Gedim::MeshMatricesDAO domainMesh(domainMeshData);

    const vector<double> baseCoordinates = geometryUtilities.EquispaceCoordinates(config.MeshParameter(), true);

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBase,
                                      rectangleHeight,
                                      baseCoordinates,
                                      baseCoordinates,
                                      domainMesh);

    const unsigned int n_1D_points = baseCoordinates.size();
    const unsigned int n_1D_squares = n_1D_points - 1;
    const unsigned int n_2D_points = domainMesh.Cell0DTotalNumber();
    const unsigned int n_2D_squares = domainMesh.Cell2DTotalNumber();

    PrintMessage(rank,
                 cerr,
                 "n_1D_points: " +
                 to_string(n_1D_points) + " " +
                 "n_1D_squares: " +
                 to_string(n_1D_squares) + " " +
                 "n_2D_points: " +
                 to_string(n_2D_points) + " " +
                 "n_2D_squares: " +
                 to_string(n_2D_squares)
                 , true);

    if (n_1D_squares % 2 > 0)
    {
      PrintMessage(rank, cerr, "n_1D_squares shall be even", true);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    PrintMessage(rank, cout, "Create Domain Mesh SUCCESS", false);

    // Export the domain mesh
    PrintMessage(rank, cout, "Export Domain Mesh...", true);
    if (rank == 0)
    {
      meshUtilities.ExportMeshToVTU(domainMesh,
                                    exportVtuFolder,
                                    "Domain_Mesh");
    }
    PrintMessage(rank, cout, "Export Domain Mesh SUCCESS", true);

    PrintMessage(rank, cout, "Compute domain geometric properties...", false);
    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        domainMesh);

    PrintMessage(rank, cout, "Compute domain geometric properties SUCCESS", false);

    /// Split mesh among processes
    const unsigned int n_1D_domains = sqrt(n_domains);

    PrintMessage(rank,
                 cerr,
                 "n_domains: " +
                 to_string(n_domains) + " " +
                 "n_1D_domains: " +
                 to_string(n_1D_domains)
                 , true);

    const unsigned int n_1D_squares_domain = n_1D_squares / n_1D_domains;
    const unsigned int n_1D_points_domain = n_1D_squares_domain + 1;

    PrintMessage(rank,
                 cerr,
                 "n_1D_points_domain: " +
                 to_string(n_1D_points_domain) + " " +
                 "n_1D_squares_domain: " +
                 to_string(n_1D_squares_domain)
                 , true);

    // Export the local domain mesh
    PrintMessage(rank, cout, "Export Local Domain Mesh...", true);

    DD_Utilities::ExportDomainToVtu(rank,
                                    n_1D_points,
                                    n_1D_squares,
                                    n_1D_domains,
                                    n_1D_squares_domain,
                                    domainMesh,
                                    exportVtuDomainFolder);
    PrintMessage(rank, cout, "Export Local Domain Mesh SUCCESS", true);

    /// Assemble System
    PrintMessage(rank, cout, "Assemble System FEM...", false);

    double numDofs = 0;

    Gedim::Eigen_SparseArray<> globalMatrixA;
    Gedim::Eigen_Array<> rightHandSide;
    Gedim::Eigen_Array<> solution;

    if (numDofs > 0)
    {
      globalMatrixA.SetSize(numDofs, numDofs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
      rightHandSide.SetSize(numDofs);
      solution.SetSize(numDofs);
    }

    PrintMessage(rank, cout, "Assemble System FEM SUCCESS", false);

    /// Solve
    PrintMessage(rank, cout, "Solve...", false);

    if (numDofs > 0)
    {
      Gedim::Eigen_CholeskySolver<> choleskySolver;
      choleskySolver.Initialize(globalMatrixA,
                                rightHandSide,
                                solution);
      choleskySolver.Solve();
    }

    PrintMessage(rank, cout, "Solve SUCCESS", false);
  }
  // ***************************************************************************
}
