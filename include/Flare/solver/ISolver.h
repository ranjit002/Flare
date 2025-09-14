#pragma once
#include "../boundary/IBoundary.h"
#include "../fluid/Fluid.h"

namespace solver
{

/**
 * @brief Abstract interface for fluid solvers.
 *
 * Solver operates on a `fluid::Fluid` object, advancing its state forward in
 * time.
 *
 * ### Responsibilities:
 * - Implementing the core simulation step through step().
 * - Managing boundary conditions (BCs) (e.g., walls, inflows, obstacles).
 * - Updates the velocity, density, and other scalar fields in the simulation.
 *
 * @note
 * - This class manages a list of BC objects, each of which
 *   defines how the fluid interacts with walls, obstacles, or inflows.
 * - Boundary objects are stored using std::unique_ptr, ensuring that
 *   `ISolver` owns their lifetime.
 */
class ISolver
{
   public:
    virtual ~ISolver() = default;

    /**
     * @brief Advances the simulation by a single timestep.
     *
     * @param fluid Fluid container object (see Flare/fluid/Fluid.h)
     * @param dt The timestep size in seconds.
     */
    virtual void step(fluid::Fluid& fluid, float dt) = 0;

    /**
     * @brief Adds new BC object to solver.
     *
     * Function transfers ownership of BC object
     * to the solver. Solver will apply all registered BCs
     * during each call to step().
     *
     * @param bc Unique pointer to BC object (IBoundary)
     *
     * @note
     * - Ownership is transferred, DON'G use the pointer after
     * passing it!
     * - Can add multiple BCs objectt.
     */
    void addBC(std::unique_ptr<boundary::IBoundary> bc)
    {
        bcs_.push_back(std::move(bc));
    }

    void addBCs(
        std::initializer_list<std::unique_ptr<boundary::IBoundary>> list)
    {
        for (auto& bc : list)
        {
            bcs_.push_back(std::move(
                const_cast<std::unique_ptr<boundary::IBoundary>&>(bc)));
        }
    }

    /**
     * @brief Read-only access to registered BCs.
     *
     * @return Const reference to vector of BCs
     */
    const std::vector<std::unique_ptr<boundary::IBoundary>>& boundaries() const
    {
        return bcs_;
    }

   protected:
    std::vector<std::unique_ptr<boundary::IBoundary>> bcs_;
};

}  // namespace solver
