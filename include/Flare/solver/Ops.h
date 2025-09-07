#pragma once
#include <Eigen/Core>

#include "Flare/fluid/IFluid.h"

namespace solver::ops
{

void addForces(fluid::IFluid& fluid, const Eigen::Vector3f& f, float dt);
void advect(fluid::IFluid& fluid, float dt);
void diffuse(fluid::IFluid& fluid, float dt, float visc = 0.1f, int iter = 10);
void project(fluid::IFluid& fluid, int iter = 20);
void advectDensity(fluid::IFluid& fluid, float dt);
void diffuseDensity(fluid::IFluid& fluid,
    float dt,
    float diff = 0.001f,
    int iter = 10);

}  // namespace solver::ops
