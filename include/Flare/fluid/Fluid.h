#pragma once
#include <Eigen/Dense>
#include <iostream>

#include "Field.h"

namespace fluid
{

class Fluid
{
 public:
  explicit Fluid(int width, int height, int depth)
      : width_(width),
        height_(height),
        depth_(depth),
        density_(width, height, depth, 0.0f),
        velocity_(width, height, depth, Eigen::Vector3f::Zero())
  {
  }

  [[nodiscard]] int width() const noexcept { return width_; }
  [[nodiscard]] int height() const noexcept { return height_; }
  [[nodiscard]] int depth() const noexcept { return depth_; }
  [[nodiscard]] size_t size() const noexcept { return density_.size(); }

  [[nodiscard]] float getDensity(int x, int y, int z) const
  {
    return density_(x, y, z);
  }

  void setDensity(int x, int y, int z, float value)
  {
    density_.set(x, y, z, value);
  }

  [[nodiscard]] const FieldFloat& density() const noexcept { return density_; }
  [[nodiscard]] FieldFloat& density() noexcept { return density_; }

  [[nodiscard]] const Eigen::Vector3f& getVelocity(int x, int y, int z) const
  {
    return velocity_(x, y, z);
  }

  void setVelocity(int x, int y, int z, const Eigen::Vector3f& vel)
  {
    velocity_.set(x, y, z, vel);
  }

  [[nodiscard]] const FieldVector& velocity() const noexcept
  {
    return velocity_;
  }
  [[nodiscard]] FieldVector& velocity() noexcept { return velocity_; }

  void clearDensity(float value = 0.0f) { density_.fill(value); }

  void clearVelocity(const Eigen::Vector3f& value = Eigen::Vector3f::Zero())
  {
    velocity_.fill(value);
  }

  void printDensity() { std::cout << density_.get().transpose() << std::endl; }

  void reset()
  {
    clearDensity();
    clearVelocity();
  }

 private:
  int width_, height_, depth_;
  FieldFloat density_;
  FieldVector velocity_;
};

}  // namespace fluid
