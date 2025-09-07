add_executable(demo_example main.cpp)

target_link_libraries(demo_example PRIVATE Flare)
target_include_directories(demo_example PRIVATE ${CMAKE_SOURCE_DIR}/include)
