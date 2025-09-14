#include "Flare/boundary/Box.h"

namespace boundary
{

bool boundary::BoxBoundary::isSolid(int x, int y, int z) const
{
    return x == 0 || x == width_ - 1 || y == 0 || y == height_ - 1 || z == 0 ||
           z == depth_ - 1;
}

Eigen::Vector3f BoxBoundary::wallVelocity(int /*x*/, int /*y*/, int /*z*/) const
{
    return Eigen::Vector3f::Zero();
}

}  // namespace boundary