#pragma once
#include <Eigen/Core>

#include "../boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"

namespace solver::ops
{

void addForces(fluid::Fluid& fluid,
    const Eigen::Vector3f& f,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);
void advect(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);
void diffuse(fluid::Fluid& fluid,
    float dt,
    float visc,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);
void project(fluid::Fluid& fluid,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);
void advectDensity(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);
void diffuseDensity(fluid::Fluid& fluid,
    float dt,
    float diff,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

}  // namespace solver::ops
