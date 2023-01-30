#include "DD_Utilities.hpp"

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
                                                                const unsigned int& n_1D_squares) const
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
}
