#ifndef __EllipticProblem_H
#define __EllipticProblem_H

#include "Configurations.hpp"

#include <string>

namespace METODI_FEM_2D
{
  class EllipticProblem_ProgramConfiguration final
  {
    public:
      EllipticProblem_ProgramConfiguration();

      inline std::string ExportFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }
      inline double GeometricTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }
      inline double MeshMaximumTriangleArea() const
      { return Gedim::Configurations::GetPropertyValue<double>("MeshMaximumTriangleArea"); }
  };

  class EllipticProblem final
  {
    private:
      const EllipticProblem_ProgramConfiguration& config;

    public:
      EllipticProblem(const EllipticProblem_ProgramConfiguration& config);
      ~EllipticProblem();

      void Run();
  };

}

#endif
