#include "Flare/solver/BasicSolver.h"

#include "Flare/solver/Ops.h"

namespace solver
{

void BasicSolver::step(fluid::Fluid& fluid, float dt)
{
    using namespace solver::ops;
    addForces(fluid, Eigen::Vector3f::Zero(), dt, bcs_);
    advect(fluid, dt, bcs_);
    diffuse(fluid, dt, viscosity_, diffuse_iter_, bcs_);
    project(fluid, project_iter_, bcs_);

    advectDensity(fluid, dt, bcs_);
    diffuseDensity(fluid, dt, diffusion_, diffuse_iter_, bcs_);
}

}  // namespace solver