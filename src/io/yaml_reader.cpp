#include "io/yaml_reader.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace {

void RequireNode(const YAML::Node& node, const std::string& name) {
  if (!node) {
    throw std::runtime_error("Missing required YAML node: " + name);
  }
}

template <typename T>
T ReadScalar(const YAML::Node& node, const std::string& name) {
  RequireNode(node, name);
  try {
    return node.as<T>();
  } catch (const std::exception& error) {
    throw std::runtime_error("Failed to parse YAML field '" + name +
                             "': " + error.what());
  }
}

}  // namespace

ProjectConfig YamlReader::ReadProjectConfig(const std::string& file_path) {
  YAML::Node root;

  try {
    root = YAML::LoadFile(file_path);
  } catch (const std::exception& error) {
    throw std::runtime_error("Failed to open YAML file '" + file_path +
                             "': " + error.what());
  }

  ProjectConfig config;

  const YAML::Node project = root["project"];
  RequireNode(project, "project");
  config.project.name = ReadScalar<std::string>(project["name"], "project.name");
  config.project.run_id =
      ReadScalar<std::string>(project["run_id"], "project.run_id");

  const YAML::Node geometry = root["geometry"];
  RequireNode(geometry, "geometry");
  config.geometry.formulation =
      ReadScalar<std::string>(geometry["formulation"], "geometry.formulation");

  const YAML::Node vessel = geometry["vessel"];
  RequireNode(vessel, "geometry.vessel");
  config.geometry.vessel.type =
      ReadScalar<std::string>(vessel["type"], "geometry.vessel.type");
  config.geometry.vessel.radius_m =
      ReadScalar<double>(vessel["radius_m"], "geometry.vessel.radius_m");
  config.geometry.vessel.height_m =
      ReadScalar<double>(vessel["height_m"], "geometry.vessel.height_m");

  const YAML::Node reflector = geometry["reflector"];
  RequireNode(reflector, "geometry.reflector");
  config.geometry.reflector.enabled =
      ReadScalar<bool>(reflector["enabled"], "geometry.reflector.enabled");
  config.geometry.reflector.type =
      ReadScalar<std::string>(reflector["type"], "geometry.reflector.type");
  config.geometry.reflector.focal_length_m =
      ReadScalar<double>(reflector["focal_length_m"],
                         "geometry.reflector.focal_length_m");
  config.geometry.reflector.vertex_z_m =
      ReadScalar<double>(reflector["vertex_z_m"],
                         "geometry.reflector.vertex_z_m");
  config.geometry.reflector.aperture_radius_m =
      ReadScalar<double>(reflector["aperture_radius_m"],
                         "geometry.reflector.aperture_radius_m");

  const YAML::Node source_patch = geometry["source_patch"];
  RequireNode(source_patch, "geometry.source_patch");
  config.geometry.source_patch.type =
      ReadScalar<std::string>(source_patch["type"], "geometry.source_patch.type");
  config.geometry.source_patch.r_min_m =
      ReadScalar<double>(source_patch["r_min_m"], "geometry.source_patch.r_min_m");
  config.geometry.source_patch.r_max_m =
      ReadScalar<double>(source_patch["r_max_m"], "geometry.source_patch.r_max_m");
  config.geometry.source_patch.z_m =
      ReadScalar<double>(source_patch["z_m"], "geometry.source_patch.z_m");

  const YAML::Node bubble_position = geometry["bubble_position"];
  RequireNode(bubble_position, "geometry.bubble_position");
  config.geometry.bubble_position.r_m =
      ReadScalar<double>(bubble_position["r_m"], "geometry.bubble_position.r_m");
  config.geometry.bubble_position.z_m =
      ReadScalar<double>(bubble_position["z_m"], "geometry.bubble_position.z_m");

  const YAML::Node mesh = root["mesh"];
  RequireNode(mesh, "mesh");
  config.mesh.h_bulk_m = ReadScalar<double>(mesh["h_bulk_m"], "mesh.h_bulk_m");
  config.mesh.h_source_m =
      ReadScalar<double>(mesh["h_source_m"], "mesh.h_source_m");
  config.mesh.h_reflector_m =
      ReadScalar<double>(mesh["h_reflector_m"], "mesh.h_reflector_m");
  config.mesh.h_bubble_zone_m =
      ReadScalar<double>(mesh["h_bubble_zone_m"], "mesh.h_bubble_zone_m");

  const YAML::Node liquid = root["liquid"];
  RequireNode(liquid, "liquid");
  config.liquid.density_kg_m3 =
      ReadScalar<double>(liquid["density_kg_m3"], "liquid.density_kg_m3");
  config.liquid.sound_speed_m_s =
      ReadScalar<double>(liquid["sound_speed_m_s"], "liquid.sound_speed_m_s");
  config.liquid.viscosity_pa_s =
      ReadScalar<double>(liquid["viscosity_pa_s"], "liquid.viscosity_pa_s");
  config.liquid.surface_tension_n_m =
      ReadScalar<double>(liquid["surface_tension_n_m"],
                         "liquid.surface_tension_n_m");
  config.liquid.static_pressure_pa =
      ReadScalar<double>(liquid["static_pressure_pa"], "liquid.static_pressure_pa");
  config.liquid.temperature_k =
      ReadScalar<double>(liquid["temperature_k"], "liquid.temperature_k");

  const YAML::Node acoustics = root["acoustics"];
  RequireNode(acoustics, "acoustics");
  config.acoustics.formulation =
      ReadScalar<std::string>(acoustics["formulation"], "acoustics.formulation");
  config.acoustics.frequency_hz =
      ReadScalar<double>(acoustics["frequency_hz"], "acoustics.frequency_hz");
  config.acoustics.time_convention =
      ReadScalar<std::string>(acoustics["time_convention"],
                              "acoustics.time_convention");

  const YAML::Node acoustic_source = acoustics["source"];
  RequireNode(acoustic_source, "acoustics.source");
  config.acoustics.source.boundary_tag =
      ReadScalar<int>(acoustic_source["boundary_tag"],
                      "acoustics.source.boundary_tag");
  config.acoustics.source.normal_velocity_m_s =
      ReadScalar<double>(acoustic_source["normal_velocity_m_s"],
                         "acoustics.source.normal_velocity_m_s");
  config.acoustics.source.mode =
      ReadScalar<std::string>(acoustic_source["mode"], "acoustics.source.mode");

  const YAML::Node rigid_boundaries = acoustics["rigid_boundaries"];
  RequireNode(rigid_boundaries, "acoustics.rigid_boundaries");
  if (!rigid_boundaries.IsSequence()) {
    throw std::runtime_error(
        "Field 'acoustics.rigid_boundaries' must be a sequence");
  }
  config.acoustics.rigid_boundaries.clear();
  for (const auto& value : rigid_boundaries) {
    config.acoustics.rigid_boundaries.push_back(value.as<int>());
  }

  const YAML::Node output = acoustics["output"];
  RequireNode(output, "acoustics.output");
  config.acoustics.output.sample_bubble_point =
      ReadScalar<bool>(output["sample_bubble_point"],
                       "acoustics.output.sample_bubble_point");
  config.acoustics.output.compute_focus_region =
      ReadScalar<bool>(output["compute_focus_region"],
                       "acoustics.output.compute_focus_region");

  return config;
}