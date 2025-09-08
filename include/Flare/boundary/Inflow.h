#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

class InflowBoundary : public IBoundary
{
 public:
  InflowBoundary(const Eigen::Vector3f& inflowVel) : inflowVelocity(inflowVel)
  {
  }

  bool isSolid(int x, int y, int z) const override
  {
    (void)x;
    (void)y;
    (void)z;
    // Left edge acts as a boundary
    return (x == 0);
  }

  Eigen::Vector3f wallVelocity(int x, int y, int z) const override
  {
    (void)x;
    (void)y;
    (void)z;
    // Return constant inflow velocity
    return inflowVelocity;
  }

 private:
  Eigen::Vector3f inflowVelocity;
};

}  // namespace boundary
