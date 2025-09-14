#include "Flare/solver/Ops.h"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"

namespace solver::ops
{

inline bool isSolidCell(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return true;
    return false;
}

inline Eigen::Vector3f getWallVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return bc->wallVelocity(x, y, z);
    return Eigen::Vector3f::Zero();
}

inline float getWallDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
    for (auto& bc : bcs)
        if (bc->isSolid(x, y, z))
            return bc->wallDensity(x, y, z);
    return 0.f;
}

/**
 * @brief Trilinearly interpolates the velocity field at a given positionition.
 *
 * Performs semi-Lagrangian backtracing by sampling the velocity field
 * at a floating-point positionition using trilinear interpolation.
 *
 * If any of the eight surrounding cells are marked as solid, the function
 * returns the wall velocity instead to enforce boundary conditions.
 *
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 * @param position
 * @param width
 * @param height
 * @param depth
 * @param oldVel Current velocity at position
 * @param idx Interpolated velocity vector at given positionition
 * @return Eigen::Vector3f interpolated velocity vector
 */
inline Eigen::Vector3f sampleVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& position,
    int width,
    int height,
    int depth,
    const auto& oldVel,
    const auto idx)
{
    int x0 = std::clamp(int(position[0]), 0, width - 1);
    int x1 = std::clamp(x0 + 1, 0, width - 1);
    int y0 = std::clamp(int(position[1]), 0, height - 1);
    int y1 = std::clamp(y0 + 1, 0, height - 1);
    int z0 = std::clamp(int(position[2]), 0, depth - 1);
    int z1 = std::clamp(z0 + 1, 0, depth - 1);

    // If any corner is solid, return wall velocity
    if (isSolidCell(bcs, x0, y0, z0) || isSolidCell(bcs, x1, y0, z0) ||
        isSolidCell(bcs, x0, y1, z0) || isSolidCell(bcs, x1, y1, z0) ||
        isSolidCell(bcs, x0, y0, z1) || isSolidCell(bcs, x1, y0, z1) ||
        isSolidCell(bcs, x0, y1, z1) || isSolidCell(bcs, x1, y1, z1))
    {
        return getWallVelocity(bcs, x0, y0, z0);
    }

    float xd = position[0] - x0;
    float yd = position[1] - y0;
    float zd = position[2] - z0;

    Eigen::Vector3f c00 =
        oldVel(idx(x0, y0, z0)) * (1 - xd) + oldVel(idx(x1, y0, z0)) * xd;
    Eigen::Vector3f c01 =
        oldVel(idx(x0, y0, z1)) * (1 - xd) + oldVel(idx(x1, y0, z1)) * xd;
    Eigen::Vector3f c10 =
        oldVel(idx(x0, y1, z0)) * (1 - xd) + oldVel(idx(x1, y1, z0)) * xd;
    Eigen::Vector3f c11 =
        oldVel(idx(x0, y1, z1)) * (1 - xd) + oldVel(idx(x1, y1, z1)) * xd;

    Eigen::Vector3f c0 = c00 * (1 - yd) + c10 * yd;
    Eigen::Vector3f c1 = c01 * (1 - yd) + c11 * yd;

    return c0 * (1 - zd) + c1 * zd;
}

/**
 * @brief Trilinearly interpolates the density field at a given positionition.
 *
 * This function is used during semi-Lagrangian advection of density.
 * It samples the density field at a position by interpolating
 * the values of the eight nearest grid cells. If any of these surrounding cells
 * are marked as solid, their density is replaced with the corresponding wall
 * density value provided by the boundary condition objects.
 *
 * @param bcs List of boundary condition objects to enforce (see
 * Flare/boundary/IBoundary.h)
 * @param position
 * @param width
 * @param height
 * @param depth
 * @param oldDensity Current density at position
 * @param idx
 * @return float interpolated density value
 */
inline float sampleDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& position,
    int width,
    int height,
    int depth,
    const auto& oldDensity,
    const auto& idx)
{
    int x0 = std::clamp(int(position[0]), 0, width - 1);
    int x1 = std::clamp(x0 + 1, 0, width - 1);
    int y0 = std::clamp(int(position[1]), 0, height - 1);
    int y1 = std::clamp(y0 + 1, 0, height - 1);
    int z0 = std::clamp(int(position[2]), 0, depth - 1);
    int z1 = std::clamp(z0 + 1, 0, depth - 1);

    auto getVal = [&](int xi, int yi, int zi)
    {
        return isSolidCell(bcs, xi, yi, zi) ? getWallDensity(bcs, xi, yi, zi)
                                            : oldDensity[idx(xi, yi, zi)];
    };

    float xd = position[0] - x0;
    float yd = position[1] - y0;
    float zd = position[2] - z0;

    float c00 = getVal(x0, y0, z0) * (1 - xd) + getVal(x1, y0, z0) * xd;
    float c01 = getVal(x0, y0, z1) * (1 - xd) + getVal(x1, y0, z1) * xd;
    float c10 = getVal(x0, y1, z0) * (1 - xd) + getVal(x1, y1, z0) * xd;
    float c11 = getVal(x0, y1, z1) * (1 - xd) + getVal(x1, y1, z1) * xd;

    float c0 = c00 * (1 - yd) + c10 * yd;
    float c1 = c01 * (1 - yd) + c11 * yd;

    return c0 * (1 - zd) + c1 * zd;
}

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