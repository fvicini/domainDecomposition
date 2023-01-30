#ifndef __EllipticProblem_H
#define __EllipticProblem_H

#include "Configurations.hpp"

#include <string>

namespace DOMAIN_DECOMPOSITION
{
  class EllipticProblem_ProgramConfiguration final
  {
    public:
      EllipticProblem_ProgramConfiguration();

      inline std::string ExportFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }
      inline double GeometricTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }
      inline double MeshParameter() const
      { return Gedim::Configurations::GetPropertyValue<double>("MeshParameter"); }
  };

  class EllipticProblem final
  {
    private:
      const EllipticProblem_ProgramConfiguration& config;

      static void PrintMessage(const int& rank,
                               ostream& output,
                               const string& message,
                               const bool& onlyMaster);
    public:
      EllipticProblem(const EllipticProblem_ProgramConfiguration& config);
      ~EllipticProblem();

      void Run(const int& rank);
  };

}

#endif
