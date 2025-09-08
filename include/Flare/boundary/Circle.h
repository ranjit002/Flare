#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

class CircleBoundary : public IBoundary
{
 public:
  CircleBoundary(Eigen::Vector3f centre, float radius)
      : centre_{centre}, radius_{radius}
  {
  }

  bool isSolid(int x, int y, int z) const override
  {
    float dist = (centre_ - Eigen::Vector3f{static_cast<float>(x),
                                static_cast<float>(y),
                                static_cast<float>(z)})
                     .norm();
    return (dist < radius_);
  }

  Eigen::Vector3f wallVelocity(int x, int y, int z) const override
  {
    return {0, 0, 0};
  }

 private:
  Eigen::Vector3f centre_;
  float radius_;
};

}  // namespace boundary