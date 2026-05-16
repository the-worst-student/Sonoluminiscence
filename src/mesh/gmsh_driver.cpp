#include "mesh/gmsh_driver.hpp"

#include "geometry/reflector_geometry.hpp"
#include "geometry/vessel_geometry.hpp"
#include "mesh/mesh_tags.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmsh.h>

namespace {

int AddPoint(double x, double y, double mesh_size) {
    return gmsh::model::occ::addPoint(x, y, 0.0, mesh_size);
}

void EnsurePointInsideVessel(const VesselGeometry& vessel_geometry,
                             double r_m,
                             double z_m,
                             const std::string& message) {
    if (r_m < 0.0) {
        throw std::invalid_argument(message);
    }

    const double vessel_radius_at_z = vessel_geometry.RadiusAtZ(z_m);

    if (r_m > vessel_radius_at_z + 1e-10) {
        throw std::invalid_argument(message);
    }
}

}  // namespace

GmshDriver::GmshDriver(const ProjectConfig& config) : config_(config) {}

void GmshDriver::BuildAxisymmetricMesh(const std::string& output_mesh_path) const {
    gmsh::initialize();
    gmsh::model::add("axisymmetric_resonator");

    VesselGeometry vessel_geometry(config_.geometry.vessel);
    const VesselGeometryData vessel_data = vessel_geometry.BuildData();

    const std::vector<VesselWallPoint>& outer_wall = vessel_data.outer_wall;

    if (outer_wall.size() < 2) {
        throw std::runtime_error("Vessel outer wall must contain at least two points");
    }

    const double vessel_height = vessel_data.height_m;
    const double top_radius = outer_wall.back().r_m;

    const double h_bulk = config_.mesh.h_bulk_m;
    const double h_reflector = config_.mesh.h_reflector_m;
    const double h_bubble = config_.mesh.h_bubble_zone_m;

    const double bubble_r = config_.geometry.bubble_position.r_m;
    const double bubble_z = config_.geometry.bubble_position.z_m;

    const double aperture = config_.geometry.reflector.aperture_radius_m;
    const double vertex_z = config_.geometry.reflector.vertex_z_m;
    const double focal_length = config_.geometry.reflector.focal_length_m;

    const double source_r_min = config_.geometry.source_patch.r_min_m;
    const double source_r_max = config_.geometry.source_patch.r_max_m;
    const double source_z = config_.geometry.source_patch.z_m;

    if (focal_length <= 0.0) {
        throw std::invalid_argument("Reflector focal_length_m must be positive");
    }

    if (aperture <= 0.0) {
        throw std::invalid_argument("Invalid reflector aperture");
    }

    if (vertex_z <= 0.0 || vertex_z > vessel_height) {
        throw std::invalid_argument("Invalid reflector vertex z");
    }

    const double reflector_top_z =
        vertex_z + aperture * aperture / (4.0 * focal_length);

    if (reflector_top_z > vessel_height + 1e-10) {
        throw std::invalid_argument("Reflector top is outside vessel height");
    }

    if (source_r_min < 0.0 || source_r_min >= source_r_max) {
        throw std::invalid_argument("Invalid source radial interval");
    }

    if (std::abs(source_z - vessel_height) > 1e-10) {
        throw std::invalid_argument(
            "For current geometry source_z must equal vessel_height");
    }

    if (source_r_max > top_radius + 1e-10) {
        throw std::invalid_argument("Source radial interval is outside vessel top radius");
    }

    if (bubble_z < 0.0 || bubble_z > vessel_height) {
        throw std::invalid_argument("Bubble z-coordinate is outside vessel bounds");
    }

    EnsurePointInsideVessel(
        vessel_geometry,
        bubble_r,
        bubble_z,
        "Bubble position is outside vessel bounds");

    const int p0 = AddPoint(0.0, 0.0, h_bulk);

    std::vector<int> outer_wall_points;
    outer_wall_points.reserve(outer_wall.size());

    for (const VesselWallPoint& point : outer_wall) {
        outer_wall_points.push_back(AddPoint(point.r_m, point.z_m, h_bulk));
    }

    const int p_bottom_outer = outer_wall_points.front();
    const int p_top_outer = outer_wall_points.back();

    const int p_source_right = AddPoint(source_r_max, vessel_height, h_bulk);
    const int p_source_left = AddPoint(source_r_min, vessel_height, h_bulk);
    const int p3 = AddPoint(0.0, vessel_height, h_bulk);

    const int l_bottom = gmsh::model::occ::addLine(p0, p_bottom_outer);

    int l_outer = -1;

    if (outer_wall_points.size() == 2) {
        l_outer = gmsh::model::occ::addLine(p_bottom_outer, p_top_outer);
    } else {
        l_outer = gmsh::model::occ::addSpline(outer_wall_points);
    }

    const int l_top_right = gmsh::model::occ::addLine(p_top_outer, p_source_right);
    const int l_source = gmsh::model::occ::addLine(p_source_right, p_source_left);
    const int l_top_left = gmsh::model::occ::addLine(p_source_left, p3);
    const int l_axis = gmsh::model::occ::addLine(p3, p0);

    const int outer_loop = gmsh::model::occ::addCurveLoop(
        {l_bottom, l_outer, l_top_right, l_source, l_top_left, l_axis});

    const int outer_surface = gmsh::model::occ::addPlaneSurface({outer_loop});

    ReflectorGeometry reflector(config_.geometry.reflector);
    const std::vector<ReflectorPoint> profile = reflector.BuildProfile(40);

    for (const ReflectorPoint& point : profile) {
        EnsurePointInsideVessel(
            vessel_geometry,
            point.r_m,
            point.z_m,
            "Reflector profile is outside vessel bounds");
    }

    EnsurePointInsideVessel(
        vessel_geometry,
        aperture,
        0.0,
        "Reflector bottom aperture point is outside vessel bounds");

    std::vector<int> reflector_points;
    reflector_points.reserve(profile.size());

    for (const ReflectorPoint& point : profile) {
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
        {l_bottom_reflector, l_right_reflector, -parabola_curve, -l_axis_reflector});

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

    const int bubble_point = AddPoint(bubble_r, bubble_z, h_bubble);

    gmsh::model::occ::synchronize();

    const int bubble_distance_field = gmsh::model::mesh::field::add("Distance");

    gmsh::model::mesh::field::setNumbers(
        bubble_distance_field,
        "PointsList",
        std::vector<double>{static_cast<double>(bubble_point)});

    const int bubble_threshold_field = gmsh::model::mesh::field::add("Threshold");

    gmsh::model::mesh::field::setNumber(
        bubble_threshold_field, "InField", bubble_distance_field);

    gmsh::model::mesh::field::setNumber(
        bubble_threshold_field, "SizeMin", h_bubble);

    gmsh::model::mesh::field::setNumber(
        bubble_threshold_field, "SizeMax", h_bulk);

    const double bubble_refine_radius_min = 0.01;
    const double bubble_refine_radius_max = 0.03;

    gmsh::model::mesh::field::setNumber(
        bubble_threshold_field, "DistMin", bubble_refine_radius_min);

    gmsh::model::mesh::field::setNumber(
        bubble_threshold_field, "DistMax", bubble_refine_radius_max);

    gmsh::model::mesh::field::setAsBackgroundMesh(bubble_threshold_field);

    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);

    gmsh::model::mesh::generate(2);
    gmsh::write(output_mesh_path);
    gmsh::finalize();
}