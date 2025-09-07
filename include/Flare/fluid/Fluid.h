#pragma once
#include "Field.h"
#include "IFluid.h"

namespace fluid
{

class Fluid : public IFluid
{
 public:
  Fluid(int width, int height, int depth)
      : width_(width),
        height_(height),
        depth_(depth),
        density(width, height, depth),
        velocity(width, height, depth)
  {
  }

  float getDensity(int x, int y, int z) const override
  {
    return density(x, y, z);
  }
  Eigen::Vector3f getVelocity(int x, int y, int z) const override
  {
    return velocity(x, y, z);
  }

  void setDensity(int x, int y, int z, float value) override
  {
    density.set(x, y, z, value);
  }
  void setVelocity(int x, int y, int z, const Eigen::Vector3f& vel) override
  {
    velocity.set(x, y, z, vel);
  }

  int width() const override { return width_; }
  int height() const override { return height_; }
  int depth() const override { return depth_; }

 private:
  int width_, height_, depth_;
  Field1D density;
  Field3D velocity;
};

}  // namespace fluid