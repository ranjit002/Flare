#pragma once
#include <Eigen/Core>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"

namespace solver::ops
{

bool isSolidCell(const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z);

Eigen::Vector3f getWallVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z);

float getWallDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z);

}  // namespace solver::ops
