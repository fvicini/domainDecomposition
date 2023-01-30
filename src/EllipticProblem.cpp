#include "EllipticProblem.hpp"

#include "Configurations.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_CholeskySolver.hpp"

using namespace std;
using namespace Eigen;

namespace METODI_FEM_2D
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
                                       numeric_limits<double>::epsilon(),
                                       "Geometric tolerance to perform 1D operations (Default: machine epsilon)");

    // Mesh parameters
    Gedim::Configurations::AddProperty("MeshMaximumTriangleArea",
                                       0.1,
                                       "Mesh 2D maximum triangle area (Default: 0.1)");
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
  void EllipticProblem::Run()
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = config.GeometricTolerance();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    /// Create folders
    const string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const string exportCsvFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    const string logFolder = exportFolder + "/Log";

    /// Export Configuration of the following Run
    Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini",
                                       false);

    /// Set Log folder
    Gedim::Output::CreateFolder(logFolder);
    Gedim::LogFile::LogFolder = logFolder;

    /// Get problem data
    const Eigen::MatrixXd domain = geometryUtilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                  1.0);

    // Export domain
    {
      Gedim::VTKUtilities vtkUtilities;
      vtkUtilities.AddPolygon(domain);
      vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("Create Domain Mesh...", true);

    Gedim::MeshMatrices domainMeshData;
    Gedim::MeshMatricesDAO domainMesh(domainMeshData);

    meshUtilities.CreateTriangularMesh(domain,
                                       config.MeshMaximumTriangleArea(),
                                       domainMesh);

    // Export the domain mesh
    meshUtilities.ExportMeshToVTU(domainMesh,
                                  exportVtuFolder,
                                  "Domain_Mesh");

    Gedim::Output::PrintStatusProgram("Create Domain Mesh");

    Gedim::Output::PrintGenericMessage("Compute domain geometric properties...", true);

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        domainMesh);

    const double h = *max_element(std::begin(meshGeometricData.Cell2DsAreas),
                                  std::end(meshGeometricData.Cell2DsAreas));

    Gedim::Output::PrintStatusProgram("Compute domain geometric properties");

    /// Assemble System
    Gedim::Output::PrintGenericMessage("Assemble System FEM...", true);

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

    Gedim::Output::PrintStatusProgram("Assemble System");

    /// Solve
    Gedim::Output::PrintGenericMessage("Solve...", true);

    if (numDofs > 0)
    {
      Gedim::Eigen_CholeskySolver<> choleskySolver;
      choleskySolver.Initialize(globalMatrixA,
                                rightHandSide,
                                solution);
      choleskySolver.Solve();
    }

    Gedim::Output::PrintStatusProgram("Solve");
  }
  // ***************************************************************************
}
