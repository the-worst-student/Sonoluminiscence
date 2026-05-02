cmake_minimum_required(VERSION 3.16)
project(sonoluminescence_resonator_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
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
        src/io/result_writer.cpp
        src/io/vtk_writer.cpp
        src/acoustics/acoustics_problem.cpp
        src/acoustics/boundary_conditions.cpp
        src/acoustics/helmholtz_solver.cpp
        src/acoustics/field_sampler.cpp
        src/acoustics/field_derivatives.cpp
        src/acoustics/gorkov_potential.cpp
        src/acoustics/gorkov_candidates.cpp
        src/coupling/pressure_coupling.cpp
        src/coupling/bubble_excitation.cpp
)

target_include_directories(solve_acoustics PRIVATE include "${GMSH_INCLUDE_DIR}")
target_link_libraries(solve_acoustics PRIVATE yaml-cpp Eigen3::Eigen ${GMSH_LIB})

add_executable(solve_bubble
        apps/solve_bubble.cpp
        src/io/yaml_reader.cpp
        src/io/csv_writer.cpp
        src/postprocess/metrics.cpp
        src/bubble/bubble_rhs.cpp
        src/bubble/gas_model.cpp
        src/bubble/ode_solver.cpp
        src/bubble/rayleigh_plesset.cpp
        src/bubble/thermal_model.cpp
        src/coupling/pressure_coupling.cpp
)

target_include_directories(solve_bubble PRIVATE include)
target_link_libraries(solve_bubble PRIVATE yaml-cpp GSL::gsl GSL::gslcblas)

add_executable(solve_bubble_candidates
        apps/solve_bubble_candidates.cpp
        src/io/bubble_excitation_reader.cpp
        src/bubble/bubble_simulation_result.cpp
        src/bubble/keller_miksis.cpp
)

target_include_directories(solve_bubble_candidates PRIVATE include)
target_link_libraries(solve_bubble_candidates PRIVATE yaml-cpp GSL::gsl GSL::gslcblas)

add_executable(run_pipeline
        apps/run_pipeline.cpp
        src/io/yaml_reader.cpp
        src/io/bubble_vtk_writer.cpp
        src/io/csv_writer.cpp
        src/io/result_writer.cpp
        src/io/vtk_writer.cpp
        src/postprocess/metrics.cpp
        src/mesh/gmsh_driver.cpp
        src/geometry/reflector_geometry.cpp
        src/acoustics/acoustics_problem.cpp
        src/acoustics/boundary_conditions.cpp
        src/acoustics/helmholtz_solver.cpp
        src/acoustics/field_sampler.cpp
        src/coupling/pressure_coupling.cpp
        src/coupling/pipeline.cpp
        src/coupling/bubble_excitation.cpp
        src/bubble/bubble_rhs.cpp
        src/bubble/gas_model.cpp
        src/bubble/ode_solver.cpp
        src/bubble/rayleigh_plesset.cpp
        src/bubble/thermal_model.cpp
)

target_include_directories(run_pipeline PRIVATE include "${GMSH_INCLUDE_DIR}")
target_link_libraries(run_pipeline PRIVATE yaml-cpp Eigen3::Eigen GSL::gsl GSL::gslcblas ${GMSH_LIB})