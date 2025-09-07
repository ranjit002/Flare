#pragma once
#include <Eigen/Core>

#include "ISolver.h"

namespace solver
{

class BasicSolver : public ISolver
{
 public:
  void step(fluid::IFluid& fluid, float dt) override;
};

}  // namespace solver
