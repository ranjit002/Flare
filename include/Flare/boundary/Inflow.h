#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

class InflowBoundary : public IBoundary
{
   public:
    explicit InflowBoundary(Eigen::Vector3f vel, float edge_density)
        : vel_(vel), edge_density_(edge_density)
    {
    }

    bool isSolid(int x, int y, int z) const override
    {
        return x == 0;  // left edge
    }

    Eigen::Vector3f wallVelocity(int x, int y, int z) const override
    {
        return vel_;
    }

    float wallDensity(int x, int y, int z) const override
    {
        return edge_density_;
    }

   private:
    Eigen::Vector3f vel_;
    float edge_density_;
};

}  // namespace boundary
