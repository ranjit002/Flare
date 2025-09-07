#include "Flare/solver/Ops.h"

#include "Flare/fluid/Field.h"
#include "Flare/fluid/IFluid.h"
namespace solver::ops
{

void addForces(fluid::IFluid& fluid, const Eigen::Vector3f& f, float dt)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  const Eigen::Vector3f deltav = (-1.0 * f) * dt;

  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        Eigen::Vector3f v = fluid.getVelocity(x, y, z);
        v += deltav;
        fluid.setVelocity(x, y, z, v);
      }
}

void advect(fluid::IFluid& fluid, float dt)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  // Make a temporary copy of velocity
  // TODO: Should add method to directly get full Field3D object fom IField
  fluid::Field3D oldVelocity(W, H, D);

  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
        oldVelocity.set(x, y, z, fluid.getVelocity(x, y, z));

  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        Eigen::Vector3f deltaPos =
            Eigen::Vector3f(x, y, z) - dt * oldVelocity(x, y, z);

        // Clamp to domain
        deltaPos[0] = std::clamp(deltaPos[0], 0.0f, float(W - 1));
        deltaPos[1] = std::clamp(deltaPos[1], 0.0f, float(H - 1));
        deltaPos[2] = std::clamp(deltaPos[2], 0.0f, float(D - 1));

        // Trilinear interpolation
        int x0 = int(deltaPos[0]), x1 = std::min(x0 + 1, W - 1);
        int y0 = int(deltaPos[1]), y1 = std::min(y0 + 1, H - 1);
        int z0 = int(deltaPos[2]), z1 = std::min(z0 + 1, D - 1);

        float xd = deltaPos[0] - x0;
        float yd = deltaPos[1] - y0;
        float zd = deltaPos[2] - z0;

        // Interpolate velocity
        Eigen::Vector3f c00 =
            oldVelocity(x0, y0, z0) * (1 - xd) + oldVelocity(x1, y0, z0) * xd;
        Eigen::Vector3f c01 =
            oldVelocity(x0, y0, z1) * (1 - xd) + oldVelocity(x1, y0, z1) * xd;
        Eigen::Vector3f c10 =
            oldVelocity(x0, y1, z0) * (1 - xd) + oldVelocity(x1, y1, z0) * xd;
        Eigen::Vector3f c11 =
            oldVelocity(x0, y1, z1) * (1 - xd) + oldVelocity(x1, y1, z1) * xd;

        Eigen::Vector3f c0 = c00 * (1 - yd) + c10 * yd;
        Eigen::Vector3f c1 = c01 * (1 - yd) + c11 * yd;

        fluid.setVelocity(x, y, z, c0 * (1 - zd) + c1 * zd);
      }
}

void diffuse(fluid::IFluid& fluid, float dt, float visc, int iter)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  fluid::Field3D oldVelocity(W, H, D);
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
        oldVelocity.set(x, y, z, fluid.getVelocity(x, y, z));

  float a = dt * visc;  // timestep-scaled viscosity

  for (int k = 0; k < iter; ++k)
  {
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          Eigen::Vector3f v =
              (fluid.getVelocity(x - 1, y, z) + fluid.getVelocity(x + 1, y, z) +
                  fluid.getVelocity(x, y - 1, z) +
                  fluid.getVelocity(x, y + 1, z) +
                  fluid.getVelocity(x, y, z - 1) +
                  fluid.getVelocity(x, y, z + 1) + a * oldVelocity(x, y, z)) /
              (6 + a);

          fluid.setVelocity(x, y, z, v);
        }
  }
}

void project(fluid::IFluid& fluid, int iter)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  fluid::Field3D velocity(W, H, D);
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
        velocity.set(x, y, z, fluid.getVelocity(x, y, z));

  fluid::Field1D div(W, H, D);  // divergence
  fluid::Field1D p(W, H, D);    // pressure

  // Compute divergence
  for (int x = 1; x < W - 1; ++x)
    for (int y = 1; y < H - 1; ++y)
      for (int z = 1; z < D - 1; ++z)
      {
        float divVal =
            0.5f * (velocity(x + 1, y, z)[0] - velocity(x - 1, y, z)[0] +
                       velocity(x, y + 1, z)[1] - velocity(x, y - 1, z)[1] +
                       velocity(x, y, z + 1)[2] - velocity(x, y, z - 1)[2]);
        div.set(x, y, z, divVal);
        p.set(x, y, z, 0.0f);
      }

  // Solve Poisson: simple Jacobi
  for (int k = 0; k < iter; ++k)
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          float newP = (p(x - 1, y, z) + p(x + 1, y, z) + p(x, y - 1, z) +
                           p(x, y + 1, z) + p(x, y, z - 1) + p(x, y, z + 1) -
                           div(x, y, z)) /
                       6.0f;
          p.set(x, y, z, newP);
        }

  // Subtract pressure gradient
  for (int x = 1; x < W - 1; ++x)
    for (int y = 1; y < H - 1; ++y)
      for (int z = 1; z < D - 1; ++z)
      {
        Eigen::Vector3f vel = velocity(x, y, z);
        vel[0] -= 0.5f * (p(x + 1, y, z) - p(x - 1, y, z));
        vel[1] -= 0.5f * (p(x, y + 1, z) - p(x, y - 1, z));
        vel[2] -= 0.5f * (p(x, y, z + 1) - p(x, y, z - 1));
        fluid.setVelocity(x, y, z, vel);
      }
}

void advectDensity(fluid::IFluid& fluid, float dt)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  fluid::Field1D oldDensity(W, H, D);
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
        oldDensity.set(x, y, z, fluid.getDensity(x, y, z));

  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
      {
        Eigen::Vector3f deltaPos =
            Eigen::Vector3f(x, y, z) - dt * fluid.getVelocity(x, y, z);

        // Clamp to domain
        deltaPos[0] = std::clamp(deltaPos[0], 0.0f, float(W - 1));
        deltaPos[1] = std::clamp(deltaPos[1], 0.0f, float(H - 1));
        deltaPos[2] = std::clamp(deltaPos[2], 0.0f, float(D - 1));

        int x0 = int(deltaPos[0]), x1 = std::min(x0 + 1, W - 1);
        int y0 = int(deltaPos[1]), y1 = std::min(y0 + 1, H - 1);
        int z0 = int(deltaPos[2]), z1 = std::min(z0 + 1, D - 1);

        float xd = deltaPos[0] - x0;
        float yd = deltaPos[1] - y0;
        float zd = deltaPos[2] - z0;

        // Trilinear interpolation for density
        float c00 =
            oldDensity(x0, y0, z0) * (1 - xd) + oldDensity(x1, y0, z0) * xd;
        float c01 =
            oldDensity(x0, y0, z1) * (1 - xd) + oldDensity(x1, y0, z1) * xd;
        float c10 =
            oldDensity(x0, y1, z0) * (1 - xd) + oldDensity(x1, y1, z0) * xd;
        float c11 =
            oldDensity(x0, y1, z1) * (1 - xd) + oldDensity(x1, y1, z1) * xd;

        float c0 = c00 * (1 - yd) + c10 * yd;
        float c1 = c01 * (1 - yd) + c11 * yd;

        fluid.setDensity(x, y, z, c0 * (1 - zd) + c1 * zd);
      }
}

void diffuseDensity(fluid::IFluid& fluid, float dt, float diff, int iter)
{
  int W = fluid.width();
  int H = fluid.height();
  int D = fluid.depth();

  fluid::Field1D oldDensity(W, H, D);
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      for (int z = 0; z < D; ++z)
        oldDensity.set(x, y, z, fluid.getDensity(x, y, z));

  float a = dt * diff;

  for (int k = 0; k < iter; ++k)
    for (int x = 1; x < W - 1; ++x)
      for (int y = 1; y < H - 1; ++y)
        for (int z = 1; z < D - 1; ++z)
        {
          float d =
              (fluid.getDensity(x - 1, y, z) + fluid.getDensity(x + 1, y, z) +
                  fluid.getDensity(x, y - 1, z) +
                  fluid.getDensity(x, y + 1, z) +
                  fluid.getDensity(x, y, z - 1) +
                  fluid.getDensity(x, y, z + 1) + a * oldDensity(x, y, z)) /
              (6 + a);

          fluid.setDensity(x, y, z, d);
        }
}

}  // namespace solver::ops