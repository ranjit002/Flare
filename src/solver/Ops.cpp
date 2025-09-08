#include "Flare/solver/Ops.h"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <vector>

#include "Flare/boundary/IBoundary.h"
#include "Flare/fluid/Fluid.h"

namespace solver::ops
{

// Check if a grid cell is solid
inline bool isSolidCell(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
  for (auto& bc : bcs)
    if (bc->isSolid(x, y, z)) return true;
  return false;
}

// Get velocity of solid cell
inline Eigen::Vector3f getWallVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
  for (auto& bc : bcs)
    if (bc->isSolid(x, y, z)) return bc->wallVelocity(x, y, z);
  return Eigen::Vector3f::Zero();
}

inline float getWallDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    int x,
    int y,
    int z)
{
  for (auto& bc : bcs)
    if (bc->isSolid(x, y, z)) return bc->wallDensity(x, y, z);
  return 0.f;
}

// Floating-point safe trilinear sampling
// Need to worry about nearby walls
inline Eigen::Vector3f sampleVelocity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& pos,
    int W,
    int H,
    int D,
    const auto& oldVel,  // can be Eigen::Matrix<Eigen::Vector3f, Dynamic, 1>
    const auto& idx)
{
  int x0 = std::clamp(int(pos[0]), 0, W - 1);
  int x1 = std::clamp(x0 + 1, 0, W - 1);
  int y0 = std::clamp(int(pos[1]), 0, H - 1);
  int y1 = std::clamp(y0 + 1, 0, H - 1);
  int z0 = std::clamp(int(pos[2]), 0, D - 1);
  int z1 = std::clamp(z0 + 1, 0, D - 1);

  // If any corner is solid, return wall velocity
  if (isSolidCell(bcs, x0, y0, z0) || isSolidCell(bcs, x1, y0, z0) ||
      isSolidCell(bcs, x0, y1, z0) || isSolidCell(bcs, x1, y1, z0) ||
      isSolidCell(bcs, x0, y0, z1) || isSolidCell(bcs, x1, y0, z1) ||
      isSolidCell(bcs, x0, y1, z1) || isSolidCell(bcs, x1, y1, z1))
  {
    return getWallVelocity(bcs, x0, y0, z0);
  }

  float xd = pos[0] - x0;
  float yd = pos[1] - y0;
  float zd = pos[2] - z0;

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

// Need to worry about nearby walls
inline float sampleDensity(
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs,
    const Eigen::Vector3f& pos,
    int W,
    int H,
    int D,
    const auto& oldDensity,
    const auto& idx)
{
  int x0 = std::clamp(int(pos[0]), 0, W - 1);
  int x1 = std::clamp(x0 + 1, 0, W - 1);
  int y0 = std::clamp(int(pos[1]), 0, H - 1);
  int y1 = std::clamp(y0 + 1, 0, H - 1);
  int z0 = std::clamp(int(pos[2]), 0, D - 1);
  int z1 = std::clamp(z0 + 1, 0, D - 1);

  auto getVal = [&](int xi, int yi, int zi)
  {
    return isSolidCell(bcs, xi, yi, zi) ? getWallDensity(bcs, xi, yi, zi)
                                        : oldDensity[idx(xi, yi, zi)];
  };

  float xd = pos[0] - x0;
  float yd = pos[1] - y0;
  float zd = pos[2] - z0;

  float c00 = getVal(x0, y0, z0) * (1 - xd) + getVal(x1, y0, z0) * xd;
  float c01 = getVal(x0, y0, z1) * (1 - xd) + getVal(x1, y0, z1) * xd;
  float c10 = getVal(x0, y1, z0) * (1 - xd) + getVal(x1, y1, z0) * xd;
  float c11 = getVal(x0, y1, z1) * (1 - xd) + getVal(x1, y1, z1) * xd;

  float c0 = c00 * (1 - yd) + c10 * yd;
  float c1 = c01 * (1 - yd) + c11 * yd;

  return c0 * (1 - zd) + c1 * zd;
}

// Add external forces
void addForces(fluid::Fluid& fluid,
    const Eigen::Vector3f& f,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  Eigen::Vector3f deltav = f * dt;
  auto& velocity = fluid.velocity().get();

  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };

#pragma omp parallel for collapse(3)
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        int i = idx(x, y, z);
        if (isSolidCell(bcs, x, y, z))
          velocity[i] = getWallVelocity(bcs, x, y, z);
        else
          velocity[i] += deltav;
      }
}

// Advect velocities
void advect(fluid::Fluid& fluid,
    float dt,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  auto oldVelocity = fluid.velocity().get();
  auto& velocity = fluid.velocity().get();

  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };

#pragma omp parallel for collapse(3)
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        Eigen::Vector3f pos =
            Eigen::Vector3f(x, y, z) - dt * oldVelocity[idx(x, y, z)];
        pos[0] = std::clamp(pos[0], 0.f, float(W - 1));
        pos[1] = std::clamp(pos[1], 0.f, float(H - 1));
        pos[2] = std::clamp(pos[2], 0.f, float(D - 1));

        velocity[idx(x, y, z)] =
            sampleVelocity(bcs, pos, W, H, D, oldVelocity, idx);
      }
}

void diffuse(fluid::Fluid& fluid,
    float dt,
    float visc,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  auto oldVelocity = fluid.velocity().get();
  auto& velocity = fluid.velocity().get();

  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };

  for (int k = 0; k < iter; ++k)
#pragma omp parallel for collapse(3)
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          if (isSolidCell(bcs, x, y, z))
          {
            velocity[idx(x, y, z)] = getWallVelocity(bcs, x, y, z);
            continue;
          }

          Eigen::Vector3f v =
              (velocity[idx(x - 1, y, z)] + velocity[idx(x + 1, y, z)] +
                  velocity[idx(x, y - 1, z)] + velocity[idx(x, y + 1, z)] +
                  velocity[idx(x, y, z - 1)] + velocity[idx(x, y, z + 1)] +
                  dt * visc * oldVelocity[idx(x, y, z)]) /
              (6 + dt * visc);

          velocity[idx(x, y, z)] = v;
        }
}

void project(fluid::Fluid& fluid,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  auto velocity = fluid.velocity().get();
  auto& liveVelocity = fluid.velocity().get();

  std::vector<float> div(W * H * D, 0.0f);
  std::vector<float> p(W * H * D, 0.0f);

  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };

// Compute divergence
#pragma omp parallel for collapse(3)
  for (int x = 1; x < W - 1; ++x)
    for (int y = 1; y < H - 1; ++y)
      for (int z = 1; z < D - 1; ++z)
      {
        int i = idx(x, y, z);
        if (isSolidCell(bcs, x, y, z)) continue;

        div[i] =
            0.5f *
            (velocity[idx(x + 1, y, z)][0] - velocity[idx(x - 1, y, z)][0] +
                velocity[idx(x, y + 1, z)][1] - velocity[idx(x, y - 1, z)][1] +
                velocity[idx(x, y, z + 1)][2] - velocity[idx(x, y, z - 1)][2]);
        p[i] = 0.0f;
      }

  // Jacobi solver for Poisson equation
  for (int k = 0; k < iter; ++k)
  {
#pragma omp parallel for collapse(3)
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          int i = idx(x, y, z);
          if (isSolidCell(bcs, x, y, z)) continue;

          p[i] = (p[idx(x - 1, y, z)] + p[idx(x + 1, y, z)] +
                     p[idx(x, y - 1, z)] + p[idx(x, y + 1, z)] +
                     p[idx(x, y, z - 1)] + p[idx(x, y, z + 1)] - div[i]) /
                 6.0f;
        }
  }
// Subtract pressure gradient
#pragma omp parallel for collapse(3)
  for (int x = 1; x < W - 1; ++x)
    for (int y = 1; y < H - 1; ++y)
      for (int z = 1; z < D - 1; ++z)
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
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  auto oldDensity = fluid.density().get();
  auto& density = fluid.density().get();
  auto& velocity = fluid.velocity().get();
  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };

#pragma omp parallel for collapse(3)
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        Eigen::Vector3f deltaPos =
            Eigen::Vector3f(x, y, z) - dt * velocity[idx(x, y, z)];
        deltaPos[0] = std::clamp(deltaPos[0], 0.f, float(W - 1));
        deltaPos[1] = std::clamp(deltaPos[1], 0.f, float(H - 1));
        deltaPos[2] = std::clamp(deltaPos[2], 0.f, float(D - 1));

        density[idx(x, y, z)] =
            sampleDensity(bcs, deltaPos, W, H, D, oldDensity, idx);
      }
}

void diffuseDensity(fluid::Fluid& fluid,
    float dt,
    float diff,
    int iter,
    const std::vector<std::unique_ptr<boundary::IBoundary>>& bcs)
{
  int W = fluid.width(), H = fluid.height(), D = fluid.depth();
  auto oldDensity = fluid.density().get();
  auto& density = fluid.density().get();
  auto idx = [&](int x, int y, int z) { return x + W * (y + H * z); };
  for (int k = 0; k < iter; ++k)
#pragma omp parallel for collapse(3)
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          int i = idx(x, y, z);
          float d = (density[idx(x - 1, y, z)] + density[idx(x + 1, y, z)] +
                        density[idx(x, y - 1, z)] + density[idx(x, y + 1, z)] +
                        density[idx(x, y, z - 1)] + density[idx(x, y, z + 1)] +
                        dt * diff * oldDensity[i]) /
                    (6 + dt * diff);

          if (isSolidCell(bcs, x, y, z))
            density[i] = getWallDensity(bcs, x, y, z);
          else
            density[i] = d;
        }
}

}  // namespace solver::ops