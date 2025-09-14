add_library(Flare STATIC    
    solver/BasicSolver.cpp    
    solver/Ops.cpp

    boundary/Box.cpp
)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type: Debug or Release" FORCE)


find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, enabling parallelism")
    target_link_libraries(Flare PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(Flare PUBLIC FLARE_USE_OPENMP)
endif()


target_include_directories(Flare
    PUBLIC
        /opt/homebrew/Cellar/eigen/3.4.0_1
        /opt/homebrew/Cellar/sfml/3.0.1/include
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)


if(MSVC)
    target_compile_options(Flare PRIVATE /W4 /permissive-)
else()
    target_compile_options(Flare PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()


target_link_directories(Flare 
    PUBLIC
        /opt/homebrew/Cellar/eigen/3.4.0_1
        /opt/homebrew/Cellar/sfml/3.0.1/lib
)


target_link_libraries(Flare
    PUBLIC Eigen3::Eigen
    PUBLIC sfml-graphics
    PUBLIC sfml-window
    PUBLIC sfml-system
)


include(GNUInstallDirs)
install(TARGETS Flare EXPORT FlareTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT FlareTargets
    FILE FlareTargets.cmake
    NAMESPACE Flare::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Flare
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
