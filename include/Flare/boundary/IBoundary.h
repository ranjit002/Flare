#pragma once
#include <Eigen/Core>

namespace boundary
{

class IBoundary
{
 public:
  virtual ~IBoundary() = default;
  virtual bool isSolid(int x, int y, int z) const = 0;
  virtual Eigen::Vector3f wallVelocity(int x, int y, int z) const = 0;
  virtual float wallDensity(int x, int y, int z) const
  {
    void(x), void(y), void(z);
    return 0.0f;
  };
};

}  // namespace boundary