#include <iostream>
#include <string>

#include "io/yaml_reader.hpp"

int main(const int argc, char** argv) {
    const std::string config_path =
        (argc > 1) ? argv[1] : "configs/base.yaml";

    try {
        const ProjectConfig config = YamlReader::ReadProjectConfig(config_path);

        std::cout << "Project: " << config.project.name << '\n';
        std::cout << "Run id: " << config.project.run_id << '\n';

        std::cout << "Vessel radius: " << config.geometry.vessel.radius_m << '\n';
        std::cout << "Vessel height: " << config.geometry.vessel.height_m << '\n';

        std::cout << "Reflector enabled: "
                  << (config.geometry.reflector.enabled ? "true" : "false") << '\n';
        std::cout << "Reflector focal length: "
                  << config.geometry.reflector.focal_length_m << '\n';

        std::cout << "Bubble position: r=" << config.geometry.bubble_position.r_m
                  << ", z=" << config.geometry.bubble_position.z_m << '\n';

        std::cout << "Frequency: " << config.acoustics.frequency_hz << '\n';
        std::cout << "Time convention: " << config.acoustics.time_convention << '\n';

        std::cout << "Liquid density: " << config.liquid.density_kg_m3 << '\n';
        std::cout << "Liquid sound speed: " << config.liquid.sound_speed_m_s << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}