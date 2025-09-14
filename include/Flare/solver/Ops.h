/**
 * @file Ops.h
 * @brief Core solver operations for the fluid simulation.
 *
 * This file contains functions for:
 * - Applying forces to the velocity field.
 * - Advection and diffusion of velocity and density.
 * - Pressure projection to enforce incompressibility.
 * - Boundary handling for solid objects.
 *
 * All functions are thread-safe and use OpenMP for parallelism.

 * @TODO: Need to eventaully split this code up to support more implementations.
 */

#pragma once
#include <Eigen/Core>

#include "../boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"

namespace solver::ops
{
/**
 * @brief Accelerate velocity vectors using a constant force boundaries are not
 * accelerated
 *
 * @note
 * - Solid cells are directly set to their wall velocity values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param force Force vector
 * @param dt Simulation timestep
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void addForces(fluid::Fluid& fluid,
    const Eigen::Vector3f& f,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

/**
 * @brief Advects the velocity field using semi-Lagrangian backtracing.
 *
 * For each cell, we trace backward in time along the velocity field
 * to find the source positionition. The velocity at this source positionition
 * is then interpolated and assigned to the current cell.
 *
 * @note
 * - Solid cells are directly set to their wall velocity values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param dt Simulation timestep
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void advectVelocity(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

/**
 * @brief Diffuses the velocity field using a simple iterative solver.
 *
 * This function simulates the diffusion of velocity over time due to viscosity.
 * It solves the discretized diffusion equation using Jacobi iteration:
 *
 * u_{new} = \frac{\sum_{neighbors} u + dt \cdot viscosity \cdot u_{old}}{6 + dt
 * \cdot viscosity}
 *
 * @note
 * - Solid cells are directly set to their wall velocity values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param dt Simulation timestep
 * @param viscosity
 * @param iter Number of Jacobi iterations to run for convergence
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void diffuseVelocity(fluid::Fluid& fluid,
    float dt,
    float visc,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

/**
 * @brief Enforces incompressibility by projecting the velocity field
 * onto a divergence-free vector field.
 *
 * The goal is to ensure the fluid is incompressible by removing its divergence
 *
 * @note
 * - Solid cells are directly set to their wall velocity values each iteration.
 * - This function uses OpenMP to parallelize over the 3D grid.
 *
 * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
 * @param iter Number of Jacobi iterations
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 */
void project(fluid::Fluid& fluid,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs);

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
