add_library(Flare STATIC    
    solver/BasicSolver.cpp    
    solver/Ops.cpp
    boundary/Wall.cpp
)

target_include_directories(Flare
    PUBLIC
        /opt/homebrew/Cellar/eigen/3.4.0_1
        /opt/homebrew/Cellar/sfml/3.0.1/include
    PUBLIC
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

# Link Dependencies
target_link_libraries(Flare
    Eigen3::Eigen
    sfml-graphics
    sfml-window
    sfml-system
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
