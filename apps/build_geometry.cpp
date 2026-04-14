#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "io/yaml_reader.hpp"
#include "mesh/gmsh_driver.hpp"

namespace {

std::string ResolveExistingPath(const std::string& path) {
    namespace fs = std::filesystem;

    const fs::path direct(path);
    if (fs::exists(direct)) {
        return direct.string();
    }

    const fs::path parent = fs::path("..") / direct;
    if (fs::exists(parent)) {
        return parent.string();
    }

    return path;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string config_path =
        argc > 1 ? ResolveExistingPath(argv[1]) : ResolveExistingPath("configs/base.yaml");
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
