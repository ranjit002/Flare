#pragma once
#include <Eigen/Core>

#include "BoundaryCondition.h"

namespace boundary
{

class Wall : public BoundaryCondition
{
 public:
  void apply(fluid::IFluid& fluid) override;
};

}  // namespace boundary