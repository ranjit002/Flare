#pragma once
#include <Eigen/Core>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"

namespace solver::ops
{

/**
 * @brief Advects the density field using semi-Lagrangian backtracing.
 *
 * For each cell, this function:
 * 1. Backtraces from the current cell center along the velocity field by `dt`.
 * 2. Finds the source position of the density using the velocity field.
 * 3. Interpolates the density at this source position using trilinear
 * interpolation.
 *
 * This method is stable even for large timesteps but introduces some numerical
 * diffusion.
 *
 * @note
 * - Solid cells are directly set to their wall density values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param dt Simulation timestep
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void advectDensity(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

/**
 * @brief Diffuses the density field using an iterative Jacobi solver.
 *
 * This function solves the diffusion equation for the density:
 *
 * d_{new} = \frac{\sum_{neighbors} d + dt \cdot diffusion \cdot d_{old}}{6 + dt
 * \cdot diff}
 *
 * @note
 * - Solid cells are directly set to their wall density values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param dt Simulation timestep
 * @param diffusion Diffusion coefficient
 * @param iter
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void diffuseDensity(fluid::Fluid& fluid,
    float dt,
    float diff,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

}  // namespace solver::ops
