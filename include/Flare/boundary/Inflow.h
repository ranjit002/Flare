#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

/**
 * @brief A boundary that injects fluid from the 'left edge' into the domain.
 *
 * Models a source at the left edge that continuously adds
 * velocity and density to the fluid. Useful for simulating
 * inflow jets or smoke sources.
 *
 */
class InflowBoundary : public IBoundary
{
   public:
    explicit InflowBoundary(Eigen::Vector3f wallVelocity, float edgeDensity)
        : wallVelocity_(wallVelocity), edgeDensity_(edgeDensity)
    {
    }

    bool isSolid(int x, int y, int z) const override
    {
        return x == 0;  // left boundary
    }

    Eigen::Vector3f wallVelocity(int x, int y, int z) const override
    {
        return wallVelocity_;
    }

    float wallDensity(int x, int y, int z) const override
    {
        return edgeDensity_;
    }

   private:
    Eigen::Vector3f wallVelocity_;
    float edgeDensity_;
};

}  // namespace boundary
