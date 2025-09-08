#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

/**
 * @brief Spherical obstacle inside fluid.
 *
 * Represents a solid circular (spherical in 3D) object.
 * Any cell whose center lies inside the radius is considered solid.
 *
 */
class CircleBoundary : public IBoundary
{
   public:
    CircleBoundary(Eigen::Vector3f centre, float radius)
        : centre_(centre), radius_(radius)
    {
    }

    bool isSolid(int x, int y, int z) const override
    {
        float dist =
            (centre_ - Eigen::Vector3f(float(x), float(y), float(z))).norm();
        return dist < radius_;
    }

    Eigen::Vector3f wallVelocity(int x, int y, int z) const override
    {
        return Eigen::Vector3f::Zero();
    }

   private:
    Eigen::Vector3f centre_;
    float radius_;
};

}  // namespace boundary