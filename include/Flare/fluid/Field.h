#pragma once
#include <Eigen/Dense>
#include <cassert>
#include <stdexcept>

namespace fluid
{
template <typename T>
class Field
{
 public:
  using StorageType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  explicit Field(int w, int h, int d, const T& default_value = T{})
      : width_(w),
        height_(h),
        depth_(d),
        data_(static_cast<Eigen::Index>(w) * h * d)
  {
    if (w <= 0 || h <= 0 || d <= 0)
    {
      throw std::invalid_argument("Field dimensions must be positive.");
    }
    data_.setConstant(default_value);
  }

  [[nodiscard]] int width() const noexcept { return width_; }
  [[nodiscard]] int height() const noexcept { return height_; }
  [[nodiscard]] int depth() const noexcept { return depth_; }
  [[nodiscard]] Eigen::Index size() const noexcept { return data_.size(); }

  void set(int x, int y, int z, const T& value)
  {
    assert(is_valid_index(x, y, z));
    data_[idx(x, y, z)] = value;
  }

  void set(const StorageType& new_data)
  {
    if (new_data.size() != data_.size())
    {
      throw std::invalid_argument("New data must match field dimensions.");
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

  [[nodiscard]] StorageType& get() noexcept { return data_; }
  [[nodiscard]] const StorageType& get() const noexcept { return data_; }

  void fill(const T& value) { data_.setConstant(value); }

 private:
  int width_, height_, depth_;
  StorageType data_;

  [[nodiscard]] inline int idx(int x, int y, int z) const noexcept
  {
    return x + width_ * (y + height_ * z);
  }

  [[nodiscard]] bool is_valid_index(int x, int y, int z) const noexcept
  {
    return (x >= 0 && x < width_) && (y >= 0 && y < height_) &&
           (z >= 0 && z < depth_);
  }
};

using FieldFloat = Field<float>;
using FieldBool = Field<bool>;
using FieldVector = Field<Eigen::Vector3f>;

}  // namespace fluid