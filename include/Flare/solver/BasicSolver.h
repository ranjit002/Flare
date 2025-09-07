#pragma once
#include <Eigen/Core>

#include "ISolver.h"

namespace solver
{

class BasicSolver : public ISolver
{
 public:
  void step(fluid::IFluid& fluid, float dt) override;
  
  float visc = 0.1f;
  float diff = 0.001f;
  int diffuse_iter = 10;
  int project_iter = 10;
};

}  // namespace solver
