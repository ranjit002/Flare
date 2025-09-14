#include "Flare/fluid/Fluid.h"

#include <iostream>

namespace fluid
{

size_t Fluid::size() const noexcept { return density_.size(); }

float Fluid::getDensity(int x, int y, int z) const { return density_(x, y, z); }

void Fluid::setDensity(int x, int y, int z, float value)
{
    density_.set(x, y, z, value);
}

const Eigen::Vector3f& Fluid::getVelocity(int x, int y, int z) const
{
    return velocity_(x, y, z);
}

void Fluid::setVelocity(int x, int y, int z, const Eigen::Vector3f& vel)
{
    velocity_.set(x, y, z, vel);
}

void Fluid::clearDensity(float value) { density_.fill(value); }

void Fluid::clearVelocity(const Eigen::Vector3f& value)
{
    velocity_.fill(value);
}

void Fluid::printDensity()
{
    std::cout << density_.get().transpose() << std::endl;
}

}  // namespace fluid