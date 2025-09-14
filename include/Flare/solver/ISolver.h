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
     * @brief Adds a single boundary condition to the solver.
     */
    template <typename BC>
    void addBC(BC&& bc)
    {
        static_assert(std::is_base_of_v<boundary::IBoundary, std::decay_t<BC>>,
            "BC must derive from boundary::IBoundary");
        bcs_.push_back(
            std::make_unique<std::decay_t<BC>>(std::forward<BC>(bc)));
    }

    /**
     * @brief Adds multiple boundary conditions in a single call.
     */
    template <typename... BCs>
    void addBCs(BCs&&... bcs)
    {
        (addBC(std::forward<BCs>(bcs)), ...);
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
