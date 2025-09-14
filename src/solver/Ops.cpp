#include "Flare/solver/Ops.h"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"
#include "Flare/solver/ops/BoundaryUtils.h"
#include "Flare/solver/ops/Samplers.h"

namespace solver::ops
{

void addForces(fluid::Fluid& fluid,
    const Eigen::Vector3f& force,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    Eigen::Vector3f deltav = force * dt;
    auto& velocity = fluid.velocity().get();

    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };

#pragma omp parallel for collapse(3)
    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            for (int z = 0; z < depth; ++z)
            {
                int i = idx(x, y, z);
                if (isSolidCell(bcs, x, y, z))
                    velocity[i] = getWallVelocity(bcs, x, y, z);
                else
                    velocity[i] += deltav;
            }
}

void advectVelocity(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    auto oldVelocity = fluid.velocity().copy();
    auto& velocity = fluid.velocity().get();

    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };

#pragma omp parallel for collapse(3)
    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            for (int z = 0; z < depth; ++z)
            {
                Eigen::Vector3f position =
                    Eigen::Vector3f(x, y, z) - dt * oldVelocity[idx(x, y, z)];
                position[0] = std::clamp(position[0], 0.f, float(width - 1));
                position[1] = std::clamp(position[1], 0.f, float(height - 1));
                position[2] = std::clamp(position[2], 0.f, float(depth - 1));

                velocity[idx(x, y, z)] = sampleVelocity(
                    bcs, position, width, height, depth, oldVelocity, idx);
            }
}

void diffuseVelocity(fluid::Fluid& fluid,
    float dt,
    float viscosity,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    auto oldVelocity = fluid.velocity().copy();
    auto& velocity = fluid.velocity().get();

    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };

    for (int k = 0; k < iter; ++k)
#pragma omp parallel for collapse(3)
        for (int x = 1; x < width - 1; ++x)
            for (int y = 1; y < height - 1; ++y)
                for (int z = 1; z < depth - 1; ++z)
                {
                    if (isSolidCell(bcs, x, y, z))
                    {
                        velocity[idx(x, y, z)] = getWallVelocity(bcs, x, y, z);
                        continue;
                    }

                    Eigen::Vector3f v =
                        (velocity[idx(x - 1, y, z)] +
                            velocity[idx(x + 1, y, z)] +
                            velocity[idx(x, y - 1, z)] +
                            velocity[idx(x, y + 1, z)] +
                            velocity[idx(x, y, z - 1)] +
                            velocity[idx(x, y, z + 1)] +
                            dt * viscosity * oldVelocity[idx(x, y, z)]) /
                        (6 + dt * viscosity);

                    velocity[idx(x, y, z)] = v;
                }
}

void project(fluid::Fluid& fluid,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
    int width = fluid.width(), height = fluid.height(), depth = fluid.depth();
    auto velocity = fluid.velocity().get();
    auto& liveVelocity = fluid.velocity().get();

    std::vector<float> div(width * height * depth, 0.0f);
    std::vector<float> p(width * height * depth, 0.0f);

    auto idx = [&](int x, int y, int z)
    { return x + width * (y + height * z); };

// Compute divergence
#pragma omp parallel for collapse(3)
    for (int x = 1; x < width - 1; ++x)
        for (int y = 1; y < height - 1; ++y)
            for (int z = 1; z < depth - 1; ++z)
            {
                int i = idx(x, y, z);
                if (isSolidCell(bcs, x, y, z))
                    continue;

                div[i] = 0.5f * (velocity[idx(x + 1, y, z)][0] -
                                    velocity[idx(x - 1, y, z)][0] +
                                    velocity[idx(x, y + 1, z)][1] -
                                    velocity[idx(x, y - 1, z)][1] +
                                    velocity[idx(x, y, z + 1)][2] -
                                    velocity[idx(x, y, z - 1)][2]);
                p[i] = 0.0f;
            }

    // Jacobi solver for Poisson equation
    for (int k = 0; k < iter; ++k)
    {
#pragma omp parallel for collapse(3)
        for (int x = 1; x < width - 1; ++x)
            for (int y = 1; y < height - 1; ++y)
                for (int z = 1; z < depth - 1; ++z)
                {
                    int i = idx(x, y, z);
                    if (isSolidCell(bcs, x, y, z))
                        continue;

                    p[i] = (p[idx(x - 1, y, z)] + p[idx(x + 1, y, z)] +
                               p[idx(x, y - 1, z)] + p[idx(x, y + 1, z)] +
                               p[idx(x, y, z - 1)] + p[idx(x, y, z + 1)] -
                               div[i]) /
                           6.0f;
                }
    }
// Subtract pressure gradient
#pragma omp parallel for collapse(3)
    for (int x = 1; x < width - 1; ++x)
        for (int y = 1; y < height - 1; ++y)
            for (int z = 1; z < depth - 1; ++z)
            {
                int i = idx(x, y, z);
                if (isSolidCell(bcs, x, y, z))
                {
                    liveVelocity[i] = getWallVelocity(bcs, x, y, z);
                    continue;
                }

                Eigen::Vector3f vel = velocity[i];
                vel[0] -= 0.5f * (p[idx(x + 1, y, z)] - p[idx(x - 1, y, z)]);
                vel[1] -= 0.5f * (p[idx(x, y + 1, z)] - p[idx(x, y - 1, z)]);
                vel[2] -= 0.5f * (p[idx(x, y, z + 1)] - p[idx(x, y, z - 1)]);
                liveVelocity[i] = vel;
            }
}

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