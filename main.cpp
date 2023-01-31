#include "Configurations.hpp"
#include "IOUtilities.hpp"

#include "EllipticProblem.hpp"

#include <mpi.h>

#define TEST_MPI false

using namespace std;

// ***************************************************************************
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n_domains;
  MPI_Comm_size(MPI_COMM_WORLD, &n_domains);

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
#if !TEST_MPI
  program.Run(rank,
              n_domains);
#else
  for (unsigned int r = 0; r < 4; r++)
  {
    program.Run(r,
                4);
  }
#endif

  // Close Program
  Gedim::Configurations::Reset();

  MPI_Finalize();

  return 0;
}
// ***************************************************************************
