#include "Flare/boundary/Wall.h"

namespace boundary
{
void apply(fluid::IFluid& fluid)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  // Zero velocity at boundaries
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
    {
      fluid.setVelocity(x, y, 0, Eigen::Vector3f::Zero());      // bottom
      fluid.setVelocity(x, y, D - 1, Eigen::Vector3f::Zero());  // top
    }

  for (int x = 0; x < W; ++x)
    for (int z = 0; z < D; ++z)
    {
      fluid.setVelocity(x, 0, z, Eigen::Vector3f::Zero());      // front
      fluid.setVelocity(x, H - 1, z, Eigen::Vector3f::Zero());  // back
    }

  for (int y = 0; y < H; ++y)
    for (int z = 0; z < D; ++z)
    {
      fluid.setVelocity(0, y, z, Eigen::Vector3f::Zero());      // left
      fluid.setVelocity(W - 1, y, z, Eigen::Vector3f::Zero());  // right
    }
}
}  // namespace boundary