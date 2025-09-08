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
    return x == 0 || x == W_ - 1 || y == 0 || y == H_ - 1 || z == 0 ||
           z == D_ - 1;
  }

  Eigen::Vector3f wallVelocity(int x, int y, int z) const override
  {
    return Eigen::Vector3f::Zero();
  }

 private:
  int W_, H_, D_;
};

}  // namespace boundary
