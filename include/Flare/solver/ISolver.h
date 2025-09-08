#pragma once
#include "../boundary/IBoundary.h"
#include "../fluid/Fluid.h"

namespace solver
{

class ISolver
{
   public:
    virtual ~ISolver() = default;
    virtual void step(fluid::Fluid& fluid, float dt) = 0;
    void addBC(std::unique_ptr<boundary::IBoundary> bc)
    {
        bcs_.push_back(std::move(bc));
    }

    const std::vector<std::unique_ptr<boundary::IBoundary>>& boundaries() const
    {
        return bcs_;
    }

   protected:
    std::vector<std::unique_ptr<boundary::IBoundary>> bcs_;
};

}  // namespace solver
