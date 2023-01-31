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
    DD_Utilities::PrintMessage(rank, cout, "Create Domain...", false);

    const Eigen::Vector3d rectangleOrigin(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBase(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeight(0.0, 1.0, 0.0);

    const Eigen::MatrixXd domain = geometryUtilities.CreateParallelogram(rectangleOrigin,
                                                                         rectangleBase,
                                                                         rectangleHeight);

    DD_Utilities::PrintMessage(rank, cout, "Create Domain SUCCESS", false);

    // Export domain
    DD_Utilities::PrintMessage(rank, cout, "Export Domain...", true);
    if (rank == 0)
    {
      Gedim::VTKUtilities vtkUtilities;
      vtkUtilities.AddPolygon(domain);
      vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }
    DD_Utilities::PrintMessage(rank, cout, "Export Domain SUCCESS", true);

    /// Create domain mesh
    DD_Utilities::PrintMessage(rank, cout,
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

    // Compute problem info
    DD_Utilities::Problem_Info problem_info = DD_Utilities::ComputeProblemInfo(rank,
                                                                               n_domains,
                                                                               domainMesh);

    DD_Utilities::PrintMessage(rank,
                               cerr,
                               "n_1D_points: " +
                               to_string(problem_info.Num_1D_points) + " " +
                               "n_1D_squares: " +
                               to_string(problem_info.Num_1D_squares) + " " +
                               "n_2D_points: " +
                               to_string(problem_info.Num_2D_points) + " " +
                               "n_2D_squares: " +
                               to_string(problem_info.Num_2D_squares)
                               , true);

    DD_Utilities::PrintMessage(rank, cout, "Create Domain Mesh SUCCESS", false);

    // Export the domain mesh
    DD_Utilities::PrintMessage(rank, cout, "Export Domain Mesh...", true);
    if (rank == 0)
    {
      meshUtilities.ExportMeshToVTU(domainMesh,
                                    exportVtuFolder,
                                    "Domain_Mesh");
    }
    DD_Utilities::PrintMessage(rank, cout, "Export Domain Mesh SUCCESS", true);

    DD_Utilities::PrintMessage(rank, cout, "Compute domain geometric properties...", false);
    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        domainMesh);

    DD_Utilities::PrintMessage(rank, cout, "Compute domain geometric properties SUCCESS", false);

    DD_Utilities::PrintMessage(rank,
                               cerr,
                               "n_domains: " +
                               to_string(n_domains) + " " +
                               "n_1D_domains: " +
                               to_string(problem_info.Num_1D_domains)
                               , true);

    DD_Utilities::PrintMessage(rank,
                               cerr,
                               "n_1D_points_domain: " +
                               to_string(problem_info.Num_1D_points_domain) + " " +
                               "n_1D_squares_domain: " +
                               to_string(problem_info.Num_1D_squares_domain)
                               , true);

    // Export the local domain mesh
    DD_Utilities::PrintMessage(rank, cout, "Export Local Domain Mesh...", false);

    DD_Utilities::ExportDomainToVtu(rank,
                                    problem_info,
                                    domainMesh,
                                    exportVtuDomainFolder);
    DD_Utilities::PrintMessage(rank, cout, "Export Local Domain Mesh SUCCESS", false);

    DD_Utilities::PrintMessage(rank, cout, "Compute DOFs...", false);
    // Compute DOFs
    DD_Utilities::DOF_Info dofs = DD_Utilities::CreateDOFs(rank,
                                                           n_domains,
                                                           problem_info,
                                                           domainMesh);
    DD_Utilities::PrintMessage(rank, cout, "Compute DOFs SUCCESS", false);

    DD_Utilities::PrintMessage(rank,
                               cerr,
                               "dofs.Num_Dirichlets: " +
                               to_string(dofs.Num_Dirichlets) + " " +
                               "dofs.Num_Internals: " +
                               to_string(dofs.Num_Internals) + " " +
                               "dofs.Num_Gamma: " +
                               to_string(dofs.Num_Gamma) + " " +
                               "dofs.Num_Globals: " +
                               to_string(dofs.Num_Globals),
                               true);

    if (rank == 0)
    {
      DD_Utilities::ExportDOFsToVtu(dofs,
                                    domainMesh,
                                    exportVtuFolder);
    }


    /// Assemble System
    DD_Utilities::PrintMessage(rank, cout, "Assemble System FEM...", false);

    const unsigned int numDofs = dofs.Num_Internals + dofs.Num_Gamma;

    Gedim::Eigen_SparseArray<> globalMatrixA;
    Gedim::Eigen_Array<> rightHandSide;
    Gedim::Eigen_Array<> internalSolution;

    if (numDofs > 0)
    {
      globalMatrixA.SetSize(numDofs, numDofs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
      rightHandSide.SetSize(numDofs);
      internalSolution.SetSize(numDofs);
    }

    DD_Utilities::Assemble(rank,
                           problem_info,
                           domainMesh,
                           meshGeometricData.Cell2DsVertices,
                           meshGeometricData.Cell2DsAreas,
                           dofs,
                           globalMatrixA,
                           rightHandSide);

    DD_Utilities::PrintMessage(rank, cout, "Assemble System FEM SUCCESS", false);

    /// Solve
    DD_Utilities::PrintMessage(rank, cout, "Solve...", false);

    if (numDofs > 0)
    {
      Gedim::Eigen_CholeskySolver<> choleskySolver;
      choleskySolver.Initialize(globalMatrixA,
                                rightHandSide,
                                internalSolution);
      choleskySolver.Solve();
    }

    DD_Utilities::PrintMessage(rank, cout, "Solve SUCCESS", false);

    DD_Utilities::PrintMessage(rank, cout, "Compute Errors...", false);

    Eigen::VectorXd errorL2_mesh, errorH1_mesh;
    DD_Utilities::ComputeErrors(rank,
                                problem_info,
                                domainMesh,
                                meshGeometricData.Cell2DsVertices,
                                meshGeometricData.Cell2DsAreas,
                                dofs,
                                internalSolution,
                                errorL2_mesh,
                                errorH1_mesh);

    const double errorL2 = sqrt(errorL2_mesh.sum());
    const double errorH1 = sqrt(errorH1_mesh.sum());

    DD_Utilities::PrintMessage(rank,
                               cerr,
                               "errorL2: " +
                               to_string(errorL2) + " " +
                               "errorH1: " +
                               to_string(errorH1),
                               true);

    // Export the local domain mesh
    DD_Utilities::PrintMessage(rank, cout, "Compute Errors SUCCESS", false);

    DD_Utilities::PrintMessage(rank, cout, "Export Local Solution...", false);

    DD_Utilities::ExportSolutionToVtu(rank,
                                      problem_info,
                                      dofs,
                                      domainMesh,
                                      internalSolution,
                                      exportVtuFolder);

    // Export the local domain mesh
    DD_Utilities::PrintMessage(rank, cout, "Export Local Solution SUCCESS", false);
  }
  // ***************************************************************************
}
