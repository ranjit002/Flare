#pragma once
#include <Eigen/Core>

namespace fluid
{

struct IFluid
{
  virtual ~IFluid() = default;

  virtual float getDensity(int x, int y, int z) const = 0;
  virtual Eigen::Vector3f getVelocity(int x, int y, int z) const = 0;

  virtual void setDensity(int x, int y, int z, float value) = 0;
  virtual void setVelocity(int x, int y, int z, const Eigen::Vector3f& vel) = 0;

  virtual int width() const = 0;
  virtual int height() const = 0;
  virtual int depth() const = 0;
};

}  // namespace fluid
