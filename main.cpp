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

  // Register program configuration
  const DOMAIN_DECOMPOSITION::EllipticProblem_ProgramConfiguration programConfiguration;

  if (rank == 0)
  {
    // Export parameters to file
    Gedim::Configurations::ExportToIni("./Parameters.ini");
  }

  // Import Parameters
  Gedim::Configurations::Initialize(argc, argv);

  DOMAIN_DECOMPOSITION::EllipticProblem program(programConfiguration);
  program.Run(rank);

  // Close Program
  Gedim::Configurations::Reset();

  MPI_Finalize();

  return 0;
}
// ***************************************************************************
