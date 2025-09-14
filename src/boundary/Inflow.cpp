#include "Flare/boundary/Inflow.h"

namespace boundary
{

bool Inflow::isSolid(int x, int /*y*/, int /*z*/) const
{
    return x == 0;  // left boundary
}

Eigen::Vector3f Inflow::wallVelocity(int /*x*/, int /*y*/, int /*z*/) const
{
    return wallVelocity_;
}

float Inflow::wallDensity(int x, int /*y*/, int /*z*/) const
{
    return x == 0 ? edgeDensity_ : 0.f;
}

}  // namespace boundary