#pragma once
#include <Eigen/Core>

#include "Flare/fluid/IFluid.h"

namespace solver::ops
{

void addForces(fluid::IFluid& fluid, const Eigen::Vector3f& f, float dt);
void advect(fluid::IFluid& fluid, float dt);
void diffuse(fluid::IFluid& fluid, float dt, float visc, int iter);
void project(fluid::IFluid& fluid, int iter);
void advectDensity(fluid::IFluid& fluid, float dt);
void diffuseDensity(fluid::IFluid& fluid, float dt, float diff, int iter);

}  // namespace solver::ops
