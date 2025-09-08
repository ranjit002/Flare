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

    bool isSolid(int x, int y, int z) const override
    {
        return x == 0 || x == width_ - 1 || y == 0 || y == height_ - 1 ||
               z == 0 || z == depth_ - 1;
    }

    Eigen::Vector3f wallVelocity(int x, int y, int z) const override
    {
        return Eigen::Vector3f::Zero();
    }

   private:
    int width_, height_, depth_;
};

}  // namespace boundary
