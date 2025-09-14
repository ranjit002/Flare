#pragma once
#include <vector>
#include <memory>
#include "Flare/boundary/IBoundary.h"
#include <Eigen/Core>

namespace solver::ops {

bool isSolidCell(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
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

} // namespace solver::ops
