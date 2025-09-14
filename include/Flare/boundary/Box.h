#pragma once
#include <Eigen/Core>

#include "IBoundary.h"

namespace boundary
{

/**
 * @brief Solid box boundary.
 *
 * Marks all outer cells at the grid edges as solid (walls) and sets velocity at
 * wall to {0, 0, 0}. Keeps fluid within simulation domain.
 *
 */
class BoxBoundary : public IBoundary
{
   public:
    BoxBoundary(int width, int height, int depth)
        : width_(width), height_(height), depth_(depth)
    {
    }

    bool isSolid(int x, int y, int z) const override;

    Eigen::Vector3f wallVelocity(int x, int y, int z) const override;

   private:
    int width_, height_, depth_;
};

}  // namespace boundary
