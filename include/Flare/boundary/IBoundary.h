#pragma once
#include <Eigen/Core>

namespace boundary
{
/**
 * @brief Interface for defining boundary conditions (BCs) in fluid simulation
 *
 * Implementations of this interface define:
 * - Which cells are considered solid (e.g., walls, obstacles).
 * - The velocity of the boundary walls (for moving or stationary boundaries).
 * - The density of walls, used when advecting or diffusing scalar quantities
 * like smoke or dye.
 *
 * Fluid solver ensures these BCs are obeyed throguhout the simulation
 *
 *
 * @note
 * - Coordinates `(x, y, z)` correspond to **cell indices**, not physical
 * positions.
 * - Derived classes must override `isSolid()` and `wallVelocity()`.
 * - `wallDensity()` has a default implementation returning 0.0f, but may be
 * overridden if you need walls to have a specific density value (e.g., to
 * simulate smoke sources).
 */
class IBoundary
{
   public:
    virtual ~IBoundary() = default;
    virtual bool isSolid(int x, int y, int z) const = 0;

    // Boundary velocity at wall
    virtual Eigen::Vector3f wallVelocity(int x, int y, int z) const = 0;
    // Boundary density at wall (defaults to 0)
    virtual float wallDensity(int x, int y, int z) const
    {
        void(x), void(y), void(z);
        return 0.0f;
    };
};

}  // namespace boundary