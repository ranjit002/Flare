#include "Flare/solver/BasicSolver.h"

#include "Flare/solver/Ops.h"

namespace solver
{

/**
 * @brief Advances the fluid state by one timestep using a basic stable fluids
 * pipeline.
 *
 * The simulation step follows the standard Navier-Stokes update sequence:
 *
 * 1. addForces - external forces.
 * 2. advect - move velocity through itself (semi-Lagrangian advection).
 * 3. diffuse - smoothing out velocity.
 * 4. project - enforce incompressibility.
 * 5. advectDensity - move scalar density through the velocity field.
 * 6. diffuseDensity - diffuse scalar density.
 *
 * Boundaries (`bcs_`) are applied at every stage.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param dt Simulation timestep
 */
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
