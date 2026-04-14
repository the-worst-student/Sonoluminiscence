cmake_minimum_required(VERSION 3.16)
project(sonoluminescence_resonator_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_library(GMSH_LIB gmsh REQUIRED)
find_path(GMSH_INCLUDE_DIR gmsh.h REQUIRED)

add_executable(build_geometry
        apps/build_geometry.cpp
        src/io/yaml_reader.cpp
        src/geometry/reflector_geometry.cpp
        src/mesh/gmsh_driver.cpp
)

target_include_directories(build_geometry PRIVATE include "${GMSH_INCLUDE_DIR}")
target_link_libraries(build_geometry PRIVATE yaml-cpp ${GMSH_LIB})

add_executable(solve_acoustics
        apps/solve_acoustics.cpp
        src/io/yaml_reader.cpp
        src/acoustics/acoustics_problem.cpp
        src/acoustics/boundary_conditions.cpp
        src/acoustics/helmholtz_solver.cpp
        src/coupling/pressure_coupling.cpp
)

target_include_directories(solve_acoustics PRIVATE include "${GMSH_INCLUDE_DIR}")
target_link_libraries(solve_acoustics PRIVATE yaml-cpp Eigen3::Eigen ${GMSH_LIB})