#include "Flare/solver/BasicSolver.h"

#include "Flare/solver/Ops.h"

namespace solver
{

void BasicSolver::step(fluid::IFluid& fluid, float dt)
{
  using namespace solver::ops;
  addForces(fluid, Eigen::Vector3f(0.f, -9.81f, 0.f), dt);
  advect(fluid, dt);
  diffuse(fluid, dt);
  project(fluid);

  advectDensity(fluid, dt);
  diffuseDensity(fluid, dt);

  applyBC(fluid);
}

}  // namespace solver