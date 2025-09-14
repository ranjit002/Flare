#pragma once
#include <Eigen/Dense>

#include "Field.h"

namespace fluid
{

/**
 * @brief Represents 3D fluid state: density + velocity.
 *
 * Provides a single access point for solver to read and modify
 * state without worrying about raw storage details.
 */
class Fluid
{
   public:
    /**
     * @brief Construct a new Fluid object
     * Initializes all density and velocity elements to the provided default
     * value.
     *
     * @param width
     * @param height
     * @param depth
     * @param default_density
     * @param default_velocity
     */
    explicit Fluid(int width,
        int height,
        int depth,
        float default_density = 0.0f,
        Eigen::Vector3f default_velocity = Eigen::Vector3f::Zero())
        : width_(width),
          height_(height),
          depth_(depth),
          density_(width, height, depth, default_density),
          velocity_(width, height, depth, default_velocity)
    {
    }

    [[nodiscard]] int width() const noexcept { return width_; }
    [[nodiscard]] int height() const noexcept { return height_; }
    [[nodiscard]] int depth() const noexcept { return depth_; }
    [[nodiscard]] size_t size() const noexcept;

    [[nodiscard]] float getDensity(int x, int y, int z) const;
    void setDensity(int x, int y, int z, float value);
    [[nodiscard]] const FieldFloat& density() const noexcept
    {
        return density_;
    }
    [[nodiscard]] FieldFloat& density() noexcept { return density_; }

    [[nodiscard]] const Eigen::Vector3f& getVelocity(int x, int y, int z) const;
    void setVelocity(int x, int y, int z, const Eigen::Vector3f& vel);
    [[nodiscard]] const FieldVector& velocity() const noexcept
    {
        return velocity_;
    }
    [[nodiscard]] FieldVector& velocity() noexcept { return velocity_; }

    void clearDensity(float value = 0.0f);
    void clearVelocity(const Eigen::Vector3f& value = Eigen::Vector3f::Zero());
    void printDensity();

   private:
    int width_, height_, depth_;
    FieldFloat density_;
    FieldVector velocity_;
};

}  // namespace fluid