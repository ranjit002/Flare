#pragma once
#include "../boundary/BoundaryCondition.h"
#include "../fluid/IFluid.h"

namespace solver
{

class ISolver
{
 public:
  virtual ~ISolver() = default;
  virtual void step(fluid::IFluid& fluid, float dt) = 0;
  void addBC(boundary::BoundaryCondition* bc) { bcs_.push_back(bc); }

  void applyBC(fluid::IFluid& fluid)
  {
    for (auto* bc : bcs_) bc->apply(fluid);
  }

 protected:
  std::vector<boundary::BoundaryCondition*> bcs_;
};

}  // namespace solver
