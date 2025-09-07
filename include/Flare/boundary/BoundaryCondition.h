#pragma once
#include "../fluid/IFluid.h"

namespace boundary
{

class BoundaryCondition
{
 public:
  virtual ~BoundaryCondition() = default;
  virtual void apply(fluid::IFluid& fluid) = 0;
};

}  // namespace boundary
