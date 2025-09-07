#pragma once
#include <Eigen/Core>

#include "ISolver.h"

namespace solver
{

class BasicSolver : public ISolver
{
 public:
  BasicSolver(float visc, float diff, int diffuse_iter, int project_iter)
      : visc_{visc},
        diff_{diff},
        diffuse_iter_{diffuse_iter},
        project_iter_{project_iter}
  {
  }
  BasicSolver()
      : visc_{0.1f}, diff_{0.001f}, diffuse_iter_{10}, project_iter_{10}
  {
  }

  void step(fluid::IFluid& fluid, float dt) override;

  float visc_;
  float diff_;
  int diffuse_iter_;
  int project_iter_;
};

}  // namespace solver
