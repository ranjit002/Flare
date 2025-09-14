#include "Flare/boundary/Circle.h"

namespace boundary
{

bool Circle::isSolid(int x, int y, int z) const
{
    float dist =
        (centre_ - Eigen::Vector3f(float(x), float(y), float(z))).norm();
    return dist < radius_;
}

Eigen::Vector3f Circle::wallVelocity(int /*x*/, int /*y*/, int /*z*/) const
{
    return Eigen::Vector3f::Zero();
}

}  // namespace boundary