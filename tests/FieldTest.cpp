#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include "Flare/fluid/Field.h"

TEST_CASE("FieldVector stores and retrieves velocity correctly",
    "[FieldVector]")
{
  const int W = 10, H = 10, D = 10;
  fluid::FieldVector field(W, H, D);

  Eigen::Vector3f velocity(1.0f, 2.0f, 3.0f);

  // Set velocity at (5,5,5)
  field.set(5, 5, 5, velocity);

  // Read it back
  Eigen::Vector3f retrieved = field(5, 5, 5);

  REQUIRE(retrieved.x() == Catch::Approx(1.0f));
  REQUIRE(retrieved.y() == Catch::Approx(2.0f));
  REQUIRE(retrieved.z() == Catch::Approx(3.0f));

  // Check default zero at another location
  Eigen::Vector3f zeroVec = field(0, 0, 0);
  REQUIRE(zeroVec.x() == Catch::Approx(0.0f));
  REQUIRE(zeroVec.y() == Catch::Approx(0.0f));
  REQUIRE(zeroVec.z() == Catch::Approx(0.0f));
}