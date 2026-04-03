#include "mesh/gmsh_driver.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmsh.h>

#include "geometry/reflector_geometry.hpp"
#include "mesh/mesh_tags.hpp"

namespace {

int AddPoint(double x, double y, double mesh_size) {
  return gmsh::model::occ::addPoint(x, y, 0.0, mesh_size);
}

}  // namespace

GmshDriver::GmshDriver(const ProjectConfig& config) : config_(config) {}

void GmshDriver::BuildAxisymmetricMesh(
    const std::string& output_mesh_path) const {
  gmsh::initialize();
  gmsh::model::add("axisymmetric_resonator");

  const double vessel_radius = config_.geometry.vessel.radius_m;
  const double vessel_height = config_.geometry.vessel.height_m;

  const double h_bulk = config_.mesh.h_bulk_m;
  const double h_reflector = config_.mesh.h_reflector_m;

  const double aperture = config_.geometry.reflector.aperture_radius_m;
  const double vertex_z = config_.geometry.reflector.vertex_z_m;
  const double focal_length = config_.geometry.reflector.focal_length_m;

  const double source_r_min = config_.geometry.source_patch.r_min_m;
  const double source_r_max = config_.geometry.source_patch.r_max_m;
  const double source_z = config_.geometry.source_patch.z_m;

  if (aperture <= 0.0 || aperture > vessel_radius) {
    throw std::invalid_argument("Invalid reflector aperture");
  }
  if (vertex_z <= 0.0 || vertex_z > vessel_height) {
    throw std::invalid_argument("Invalid reflector vertex z");
  }
  if (source_r_min < 0.0 || source_r_max > vessel_radius ||
      source_r_min >= source_r_max) {
    throw std::invalid_argument("Invalid source radial interval");
  }
  if (std::abs(source_z - vessel_height) > 1e-10) {
    throw std::invalid_argument(
        "For current geometry source_z must equal vessel_height");
  }

  const int p0 = AddPoint(0.0, 0.0, h_bulk);
  const int p1 = AddPoint(vessel_radius, 0.0, h_bulk);
  const int p2 = AddPoint(vessel_radius, vessel_height, h_bulk);
  const int p_source_right = AddPoint(source_r_max, vessel_height, h_bulk);
  const int p_source_left = AddPoint(source_r_min, vessel_height, h_bulk);
  const int p3 = AddPoint(0.0, vessel_height, h_bulk);

  const int l_bottom = gmsh::model::occ::addLine(p0, p1);
  const int l_outer = gmsh::model::occ::addLine(p1, p2);
  const int l_top_right = gmsh::model::occ::addLine(p2, p_source_right);
  const int l_source = gmsh::model::occ::addLine(p_source_right, p_source_left);
  const int l_top_left = gmsh::model::occ::addLine(p_source_left, p3);
  const int l_axis = gmsh::model::occ::addLine(p3, p0);

  const int outer_loop = gmsh::model::occ::addCurveLoop(
      {l_bottom, l_outer, l_top_right, l_source, l_top_left, l_axis});
  const int outer_surface = gmsh::model::occ::addPlaneSurface({outer_loop});

  ReflectorGeometry reflector(config_.geometry.reflector);
  const std::vector<ReflectorPoint> profile = reflector.BuildProfile(40);

  std::vector<int> reflector_points;
  reflector_points.reserve(profile.size());

  for (const auto& point : profile) {
    reflector_points.push_back(AddPoint(point.r_m, point.z_m, h_reflector));
  }

  const int parabola_curve = gmsh::model::occ::addSpline(reflector_points);

  const int p_axis_vertex = reflector_points.front();
  const int p_aperture = reflector_points.back();
  const int p_bottom_aperture = AddPoint(aperture, 0.0, h_reflector);

  const int l_axis_reflector = gmsh::model::occ::addLine(p0, p_axis_vertex);
  const int l_right_reflector =
      gmsh::model::occ::addLine(p_bottom_aperture, p_aperture);
  const int l_bottom_reflector =
      gmsh::model::occ::addLine(p0, p_bottom_aperture);

  const int reflector_loop = gmsh::model::occ::addCurveLoop(
      {l_bottom_reflector, l_right_reflector, -parabola_curve,
       -l_axis_reflector});
  const int reflector_surface =
      gmsh::model::occ::addPlaneSurface({reflector_loop});

  std::vector<std::pair<int, int>> cut_out;
  std::vector<std::vector<std::pair<int, int>>> cut_map;

  gmsh::model::occ::cut(
      {{2, outer_surface}},
      {{2, reflector_surface}},
      cut_out,
      cut_map,
      -1,
      true,
      true);

  gmsh::model::occ::synchronize();

  if (cut_out.empty()) {
    throw std::runtime_error("Failed to construct fluid surface");
  }

  int fluid_surface = -1;
  for (const auto& entity : cut_out) {
    if (entity.first == 2) {
      fluid_surface = entity.second;
      break;
    }
  }

  if (fluid_surface < 0) {
    throw std::runtime_error("Fluid surface was not found after cut");
  }

  gmsh::model::addPhysicalGroup(
      2, {fluid_surface}, static_cast<int>(SurfaceTag::cFluid));
  gmsh::model::setPhysicalName(
      2, static_cast<int>(SurfaceTag::cFluid), "fluid");

  std::vector<std::pair<int, int>> fluid_boundary;
  gmsh::model::getBoundary(
      {{2, fluid_surface}}, fluid_boundary, true, false, false);

  std::vector<int> wall_curves;
  std::vector<int> axis_curves;
  std::vector<int> source_curves;
  std::vector<int> reflector_curves;

  const double tol_axis = 1e-6;
  const double tol_top = 1e-6;
  const double tol_source = 1e-6;
  const double tol_reflector = 1e-6;

  const double reflector_top_z =
      vertex_z + aperture * aperture / (4.0 * focal_length);

  for (const auto& entity : fluid_boundary) {
    if (entity.first != 1) {
      continue;
    }

    const int tag = entity.second;

    double xmin = 0.0;
    double ymin = 0.0;
    double zmin = 0.0;
    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;

    gmsh::model::getBoundingBox(
        1, tag, xmin, ymin, zmin, xmax, ymax, zmax);

    const bool is_axis =
        std::abs(xmin) <= tol_axis && std::abs(xmax) <= tol_axis;

    const bool is_top =
        std::abs(ymin - source_z) <= tol_top &&
        std::abs(ymax - source_z) <= tol_top;

    const bool is_source =
        is_top &&
        xmin >= source_r_min - tol_source &&
        xmax <= source_r_max + tol_source;

    const bool is_reflector_parabola =
        ymin >= vertex_z - tol_reflector &&
        xmax <= aperture + tol_reflector &&
        ymax <= reflector_top_z + tol_reflector &&
        !(is_axis || is_top);

    const bool is_reflector_side =
        std::abs(xmin - aperture) <= tol_reflector &&
        std::abs(xmax - aperture) <= tol_reflector &&
        ymin >= -tol_reflector &&
        ymax <= reflector_top_z + tol_reflector;

    const bool is_reflector = is_reflector_parabola || is_reflector_side;

    if (is_source) {
      source_curves.push_back(tag);
    } else if (is_axis) {
      axis_curves.push_back(tag);
    } else if (is_reflector) {
      reflector_curves.push_back(tag);
    } else {
      wall_curves.push_back(tag);
    }
  }

  if (!wall_curves.empty()) {
    gmsh::model::addPhysicalGroup(
        1, wall_curves, static_cast<int>(BoundaryTag::cWall));
    gmsh::model::setPhysicalName(
        1, static_cast<int>(BoundaryTag::cWall), "wall");
  }

  if (!reflector_curves.empty()) {
    gmsh::model::addPhysicalGroup(
        1, reflector_curves, static_cast<int>(BoundaryTag::cReflector));
    gmsh::model::setPhysicalName(
        1, static_cast<int>(BoundaryTag::cReflector), "reflector");
  }

  if (!source_curves.empty()) {
    gmsh::model::addPhysicalGroup(
        1, source_curves, static_cast<int>(BoundaryTag::cSource));
    gmsh::model::setPhysicalName(
        1, static_cast<int>(BoundaryTag::cSource), "source");
  }

  if (!axis_curves.empty()) {
    gmsh::model::addPhysicalGroup(
        1, axis_curves, static_cast<int>(BoundaryTag::cAxis));
    gmsh::model::setPhysicalName(
        1, static_cast<int>(BoundaryTag::cAxis), "axis");
  }

  gmsh::model::mesh::generate(2);
  gmsh::write(output_mesh_path);
  gmsh::finalize();
}