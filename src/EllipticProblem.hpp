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
      inline unsigned int SchurSolverMaxIterations() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("SchurSolverMaxIterations"); }
      inline double SchurSolverTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("SchurSolverTolerance"); }
      inline bool SchurSolverType() const
      { return Gedim::Configurations::GetPropertyValue<bool>("SchurSolverType"); }
  };

  class EllipticProblem final
  {
    private:
      const EllipticProblem_ProgramConfiguration& config;

    public:
      EllipticProblem(const EllipticProblem_ProgramConfiguration& config);
      ~EllipticProblem();

      void Run(const int& rank,
               const int& n_domains) const;
  };

}

#endif
