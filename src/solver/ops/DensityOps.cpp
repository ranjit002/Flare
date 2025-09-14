#include "Flare/fluid/Fluid.h"
#include "Flare/solver/ops/BoundaryUtils.h"
#include "Flare/solver/ops/Samplers.h"

namespace solver::ops
{

void advectDensity(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    auto oldDensity = fluid.density().copy();
    auto& density = fluid.density().get();
    auto& velocity = fluid.velocity().get();
    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };

#pragma omp parallel for collapse(3)
    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            for (int z = 0; z < depth; ++z)
            {
                Eigen::Vector3f deltaPos =
                    Eigen::Vector3f(x, y, z) - dt * velocity[idx(x, y, z)];
                deltaPos[0] = std::clamp(deltaPos[0], 0.f, float(width - 1));
                deltaPos[1] = std::clamp(deltaPos[1], 0.f, float(height - 1));
                deltaPos[2] = std::clamp(deltaPos[2], 0.f, float(depth - 1));

                density[idx(x, y, z)] = sampleDensity(
                    bcs, deltaPos, width, height, depth, oldDensity, idx);
            }
}

void diffuseDensity(fluid::Fluid& fluid,
    float dt,
    float diffusion,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    auto oldDensity = fluid.density().copy();
    auto& density = fluid.density().get();
    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };
    for (int k = 0; k < iter; ++k)
#pragma omp parallel for collapse(3)
        for (int x = 1; x < width - 1; ++x)
            for (int y = 1; y < height - 1; ++y)
                for (int z = 1; z < depth - 1; ++z)
                {
                    int i = idx(x, y, z);
                    float d =
                        (density[idx(x - 1, y, z)] + density[idx(x + 1, y, z)] +
                            density[idx(x, y - 1, z)] +
                            density[idx(x, y + 1, z)] +
                            density[idx(x, y, z - 1)] +
                            density[idx(x, y, z + 1)] +
                            dt * diffusion * oldDensity[i]) /
                        (6 + dt * diffusion);

                    if (isSolidCell(bcs, x, y, z))
                        density[i] = getWallDensity(bcs, x, y, z);
                    else
                        density[i] = d;
                }
}

}  // namespace solver::ops