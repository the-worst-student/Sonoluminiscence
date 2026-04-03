#include <exception>
#include <iostream>
#include <string>

#include <io/yaml_reader.hpp>
#include <mesh/gmsh_driver.hpp>

int main(int argc, char **argv) {
    const std::string config_path = argc > 1 ? argv[1] : "config/base.yaml";
    const std::string output_mesh_path = argc > 2 ? argv[2] : "mesh.msh";

    try {
        const ProjectConfig config = YamlReader::ReadProjectConfig(config_path);

        std::cout << "Building axisymmetric mesh\n";
        std::cout << "Config: " << config_path << '\n';
        std::cout << "Output: " << output_mesh_path << '\n';

        GmshDriver driver(config);
        driver.BuildAxisymmetricMesh(output_mesh_path);

        std::cout << "Mesh successfully written to: "
                  << output_mesh_path << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error while building geometry: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}