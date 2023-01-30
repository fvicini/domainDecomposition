#include "Configurations.hpp"
#include "IOUtilities.hpp"

#include "EllipticProblem.hpp"

#include <mpi.h>

using namespace std;

// ***************************************************************************
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::cout.precision(4);
  std::cout<< std::scientific<< "Process "<< rank<< std::endl;

  //  // Register program configuration
  //  const METODI_FEM_2D::EllipticProblem_ProgramConfiguration programConfiguration;

  //  // Import Parameters
  //  if (!Gedim::Output::FileExists("./Parameters.ini"))
  //    Gedim::Configurations::ExportToIni("./Parameters.ini",
  //                                       false);
  //  else
  //    Gedim::Configurations::InitializeFromIni("./Parameters.ini");
  //  Gedim::Configurations::Initialize(argc, argv);

  //  METODI_FEM_2D::EllipticProblem program(programConfiguration);
  //  program.Run();

  //  // Close Program
  //  Gedim::Configurations::Reset();

  MPI_Finalize();

  return 0;
}
// ***************************************************************************
