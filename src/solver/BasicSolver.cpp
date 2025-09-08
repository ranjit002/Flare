#include "Flare/solver/BasicSolver.h"

#include "Flare/solver/Ops.h"

namespace solver
{

void BasicSolver::step(fluid::Fluid& fluid, float dt)
{
  using namespace solver::ops;
  addForces(fluid, Eigen::Vector3f(0.f, -9.81f, 0.f), dt, bcs_);
  advect(fluid, dt, bcs_);
  diffuse(fluid, dt, visc_, diffuse_iter_, bcs_);
  project(fluid, project_iter_, bcs_);

  advectDensity(fluid, dt);
  diffuseDensity(fluid, dt, diff_, diffuse_iter_);
}

}  // namespace solver