#pragma once
#include <Eigen/Dense>
#include <cassert>
#include <stdexcept>

namespace fluid
{

/**
 * @brief Generic 3D field class for representing fluid data.
 *
 * This class stores a 3D grid of data elements (e.g., scalars, vectors)
 * in a flat Eigen vector for cache efficiency.
 *
 * ### Memory Layout:
 * The data is stored in a flattened, row-major format:
 * index = x + width * (y + height * z)
 *
 * @tparam T The type of each field element (e.g., float, bool,
 * Eigen::Vector3f).
 */
template <typename T>
class Field
{
   public:
    /// Where data is stored
    using StorageType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    /**
     * @brief  Construct a new Field object
     *
     * Initializes all elements to the provided default value.
     *
     * @param width
     * @param height
     * @param depth
     * @param default_value The initial value for all field cells.
     *
     * @throws std::invalid_argument if any dimension is <= 0.
     */
    explicit Field(int width,
        int height,
        int depth,
        const T& default_value = T{})
        : width_(width),
          height_(height),
          depth_(depth),
          data_(static_cast<Eigen::Index>(width) * height * depth)
    {
        if (width <= 0 || height <= 0 || depth <= 0)
        {
            throw std::invalid_argument("Field dimensions must be positive.");
        }
        data_.setConstant(default_value);
    }

    [[nodiscard]] int width() const noexcept { return width_; }
    [[nodiscard]] int height() const noexcept { return height_; }
    [[nodiscard]] int depth() const noexcept { return depth_; }
    [[nodiscard]] Eigen::Index size() const noexcept { return data_.size(); }

    /**
     * @brief Sets cell to given value.
     *
     * @param x
     * @param y
     * @param z
     * @param value
     */
    void set(int x, int y, int z, const T& value)
    {
        assert(is_valid_index(x, y, z));
        data_[idx(x, y, z)] = value;
    }

    /**
     * @brief Replaces the entire field data with given field.
     *
     * @param new_data
     * @note Must match the total field size (`width * height * depth`).
     *
     * @throws std::invalid_argument if `new_data.size()` does not match.
     */
    void set(const StorageType& new_data)
    {
        if (new_data.size() != data_.size())
        {
            throw std::invalid_argument(
                "New data must match field dimensions.");
        }
        data_ = new_data;
    }

    [[nodiscard]] const T& operator()(int x, int y, int z) const
    {
        assert(is_valid_index(x, y, z));
        return data_[idx(x, y, z)];
    }

    [[nodiscard]] T& operator()(int x, int y, int z)
    {
        assert(is_valid_index(x, y, z));
        return data_[idx(x, y, z)];
    }

    /**
     * @brief Direct access to the underlying Eigen vector (mutable).
     *
     * Useful for vectorized operations.
     */
    [[nodiscard]] StorageType& get() noexcept { return data_; }

    /**
     * @brief Direct access to the underlying Eigen vector (read-only).
     */
    [[nodiscard]] const StorageType& get() const noexcept { return data_; }

    /**
     * @brief Copy of underlying Eigen vector.
     */
    [[nodiscard]] StorageType copy() const noexcept { return data_; }

    /**
     * @brief Fills the entire field with value.
     *
     * @param value
     */
    void fill(const T& value) { data_.setConstant(value); }

   private:
    int width_;
    int height_;
    int depth_;
    StorageType data_;

    /**
     * @brief Gets flat index int contigous array.
     */
    [[nodiscard]] inline int idx(int x, int y, int z) const noexcept
    {
        return x + width_ * (y + height_ * z);
    }

    /**
     * @brief Checks if a coordinate is within valid bounds.
     *
     * @param x
     * @param y
     * @param z
     * @return bool
     */
    [[nodiscard]] bool is_valid_index(int x, int y, int z) const noexcept
    {
        return (x >= 0 && x < width_) && (y >= 0 && y < height_) &&
               (z >= 0 && z < depth_);
    }
};

}  // namespace fluid
