#pragma once
#include <Eigen/Core>
#include <vector>

namespace fluid
{

class Field1D
{
 public:
  Field1D(int w, int h, int d)
      : width(w), height(h), depth(d), phi(w * h * d, 0.0f)
  {
  }

  void set(int x, int y, int z, float value) { phi[idx(x, y, z)] = value; }

  float operator()(int x, int y, int z) const { return phi[idx(x, y, z)]; }

 private:
  int width, height;
  [[maybe_unused]] int depth;
  std::vector<float> phi;

  inline int idx(int x, int y, int z) const
  {
    return x + width * (y + height * z);
  }
};

class Field3D
{
 public:
  Field3D(int w, int h, int d)
      : width(w),
        height(h),
        depth(d),
        vx(w * h * d, 0.0f),
        vy(w * h * d, 0.0f),
        vz(w * h * d, 0.0f)
  {
  }

  void set(int x, int y, int z, const Eigen::Vector3f& vector)
  {
    vx[idx(x, y, z)] = vector.x();
    vy[idx(x, y, z)] = vector.y();
    vz[idx(x, y, z)] = vector.z();
  }

  Eigen::Vector3f operator()(int x, int y, int z) const
  {
    return Eigen::Vector3f(
        vx[idx(x, y, z)], vy[idx(x, y, z)], vz[idx(x, y, z)]);
  }

 private:
  int width, height;
  [[maybe_unused]] int depth;
  std::vector<float> vx, vy, vz;

  inline int idx(int x, int y, int z) const
  {
    return x + width * (y + height * z);
  }
};

}  // namespace fluid