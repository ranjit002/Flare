#pragma once
#include <Eigen/Core>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"

namespace solver::ops
{

template <typename TVel, typename TIdx>
Eigen::Vector3f sampleVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& position,
    int width,
    int height,
    int depth,
    const TVel& oldVel,
    const TIdx& idx)
{
    int x0 = std::clamp(int(position[0]), 0, width - 1);
    int x1 = std::clamp(x0 + 1, 0, width - 1);
    int y0 = std::clamp(int(position[1]), 0, height - 1);
    int y1 = std::clamp(y0 + 1, 0, height - 1);
    int z0 = std::clamp(int(position[2]), 0, depth - 1);
    int z1 = std::clamp(z0 + 1, 0, depth - 1);

    auto isSolid = [&](int x, int y, int z)
    {
        for (auto& bc : bcs)
            if (bc->isSolid(x, y, z))
                return true;
        return false;
    };

    auto getWallVel = [&](int x, int y, int z)
    {
        for (auto& bc : bcs)
            if (bc->isSolid(x, y, z))
                return bc->wallVelocity(x, y, z);
        return Eigen::Vector3f::Zero().eval();
    };

    if (isSolid(x0, y0, z0) || isSolid(x1, y0, z0) || isSolid(x0, y1, z0) ||
        isSolid(x1, y1, z0) || isSolid(x0, y0, z1) || isSolid(x1, y0, z1) ||
        isSolid(x0, y1, z1) || isSolid(x1, y1, z1))
    {
        return getWallVel(x0, y0, z0);
    }

    float xd = position[0] - x0;
    float yd = position[1] - y0;
    float zd = position[2] - z0;

    Eigen::Vector3f c00 =
        oldVel(idx(x0, y0, z0)) * (1 - xd) + oldVel(idx(x1, y0, z0)) * xd;
    Eigen::Vector3f c01 =
        oldVel(idx(x0, y0, z1)) * (1 - xd) + oldVel(idx(x1, y0, z1)) * xd;
    Eigen::Vector3f c10 =
        oldVel(idx(x0, y1, z0)) * (1 - xd) + oldVel(idx(x1, y1, z0)) * xd;
    Eigen::Vector3f c11 =
        oldVel(idx(x0, y1, z1)) * (1 - xd) + oldVel(idx(x1, y1, z1)) * xd;

    Eigen::Vector3f c0 = c00 * (1 - yd) + c10 * yd;
    Eigen::Vector3f c1 = c01 * (1 - yd) + c11 * yd;

    return c0 * (1 - zd) + c1 * zd;
}

template <typename TDensity, typename TIdx>
float sampleDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& position,
    int width,
    int height,
    int depth,
    const TDensity& oldDensity,
    const TIdx& idx)
{
    int x0 = std::clamp(int(position[0]), 0, width - 1);
    int x1 = std::clamp(x0 + 1, 0, width - 1);
    int y0 = std::clamp(int(position[1]), 0, height - 1);
    int y1 = std::clamp(y0 + 1, 0, height - 1);
    int z0 = std::clamp(int(position[2]), 0, depth - 1);
    int z1 = std::clamp(z0 + 1, 0, depth - 1);

    auto getVal = [&](int xi, int yi, int zi)
    {
        for (auto& bc : bcs)
            if (bc->isSolid(xi, yi, zi))
                return bc->wallDensity(xi, yi, zi);
        return oldDensity[idx(xi, yi, zi)];
    };

    float xd = position[0] - x0;
    float yd = position[1] - y0;
    float zd = position[2] - z0;

    float c00 = getVal(x0, y0, z0) * (1 - xd) + getVal(x1, y0, z0) * xd;
    float c01 = getVal(x0, y0, z1) * (1 - xd) + getVal(x1, y0, z1) * xd;
    float c10 = getVal(x0, y1, z0) * (1 - xd) + getVal(x1, y1, z0) * xd;
    float c11 = getVal(x0, y1, z1) * (1 - xd) + getVal(x1, y1, z1) * xd;

    float c0 = c00 * (1 - yd) + c10 * yd;
    float c1 = c01 * (1 - yd) + c11 * yd;

    return c0 * (1 - zd) + c1 * zd;
}

}  // namespace solver::ops
