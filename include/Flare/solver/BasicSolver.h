#pragma once
#include <Eigen/Core>

#include "ISolver.h"

namespace solver
{

class BasicSolver : public ISolver
{
   public:
    BasicSolver(float viscosity, float diffusion, int diffuse_iter, int project_iter)
        : viscosity_{viscosity},
          diffusion_{diffusion},
          diffuse_iter_{diffuse_iter},
          project_iter_{project_iter}
    {
    }
    BasicSolver()
        : viscosity_{0.1f}, diffusion_{0.001f}, diffuse_iter_{10}, project_iter_{10}
    {
    }

    void step(fluid::Fluid& fluid, float dt) override;

   private:
    float viscosity_;
    float diffusion_;
    int diffuse_iter_;
    int project_iter_;
};

}  // namespace solver
