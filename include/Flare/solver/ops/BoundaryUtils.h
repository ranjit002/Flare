#pragma once
#include <Eigen/Core>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"

namespace solver::ops
{

inline bool isSolidCell(const std ::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return true;
    return false;
}

inline Eigen::Vector3f getWallVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return bc->wallVelocity(x, y, z);
    return Eigen::Vector3f::Zero();
}

inline float getWallDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return bc->wallDensity(x, y, z);
    return 0.f;
}

}  // namespace solver::ops
