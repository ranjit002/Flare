#pragma once
#include <Eigen/Core>

#include "ISolver.h"

namespace solver
{

class BasicSolver : public ISolver
{
   public:
    BasicSolver(float viscosity,
        float diffusion,
        int diffuse_iter,
        int project_iter)
        : viscosity_{viscosity},
          diffusion_{diffusion},
          diffuse_iter_{diffuse_iter},
          project_iter_{project_iter}
    {
    }
    BasicSolver()
        : viscosity_{0.1f},
          diffusion_{0.001f},
          diffuse_iter_{10},
          project_iter_{10}
    {
    }

    /**
     * @brief Advances the fluid state by one timestep using a basic stable
     * fluids pipeline.
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
    void step(fluid::Fluid& fluid, float dt) override;

   private:
    float viscosity_;
    float diffusion_;
    int diffuse_iter_;
    int project_iter_;
};

}  // namespace solver
