cmake_minimum_required(VERSION 3.20)
project(Flare VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_GENERATOR)
  set(CMAKE_GENERATOR "Ninja")
endif()

# Dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(tests)