#ifndef __DD_Utilities_H
#define __DD_Utilities_H

#include <vector>

namespace DOMAIN_DECOMPOSITION
{
  class DD_Utilities final
  {
    public:
      struct Square_Info final
      {
          unsigned int Index;
          std::vector<unsigned int> Points_Index;
      };

    private:

    public:
      DD_Utilities();
      ~DD_Utilities();

      Square_Info GetSquare_Info_Global(const unsigned int& i,
                                        const unsigned int& j,
                                        const unsigned int& n_1D_points,
                                        const unsigned int& n_1D_squares) const;
  };

}

#endif
