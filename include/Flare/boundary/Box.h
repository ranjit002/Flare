#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

class BoxBoundary : public IBoundary
{
 public:
  BoxBoundary(int W, int H, int D) : W_(W), H_(H), D_(D) {}

  bool isSolid(int x, int y, int z) const override
  {
    return x == 0 || y == 0 || z == 0 || x == W_ || y == H_ || z == D_;
  }

  Eigen::Vector3f wallVelocity(int x, int y, int z) const override
  {
    (void)x;
    (void)y;
    (void)z;
    // Return constant inflow velocity
    return {0, 0, 0};
  }

 private:
  int W_, H_, D_;
};

}  // namespace boundary
