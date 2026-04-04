#include <exception>
#include <iostream>
#include <string>

#include "acoustics/acoustics_problem.hpp"
#include "acoustics/boundary_conditions.hpp"
#include "acoustics/helmholtz_solver.hpp"
#include "io/yaml_reader.hpp"

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/base.yaml";
  const std::string mesh_path = argc > 2 ? argv[2] : "build/mesh.msh";

  try {
    const ProjectConfig config = YamlReader::ReadProjectConfig(config_path);
    const AcousticsProblem problem = AcousticsProblem::FromConfig(config);
    const BoundaryConditions conditions(problem);

    HelmholtzSolver solver(problem, conditions);
    solver.LoadMesh(mesh_path);
    solver.Assemble();
    solver.Solve();

    std::cout << "Acoustics solver finished\n";
    std::cout << "Config: " << config_path << '\n';
    std::cout << "Mesh: " << mesh_path << '\n';
    std::cout << "Mesh loaded: " << std::boolalpha << solver.IsMeshLoaded() << '\n';
    std::cout << "Assembled: " << std::boolalpha << solver.IsAssembled() << '\n';
    std::cout << "Solved: " << std::boolalpha << solver.IsSolved() << '\n';
  } catch (const std::exception& error) {
    std::cerr << "Error while solving acoustics: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
