#include <geometry/reflector_geometry.hpp>

#include <stdexcept>

ReflectorGeometry::ReflectorGeometry(const ReflectorConfig& config) : config_(config) {}

std::vector<ReflectorPoint> ReflectorGeometry::BuildProfile(int num_points) const {
    if (num_points < 2) {
        throw std::invalid_argument("Reflector profile requires at least 2 points");
    }

    std::vector<ReflectorPoint> points;
    points.reserve(num_points);

    const double aperture = config_.aperture_radius_m;
    const double focal_length = config_.focal_length_m;
    const double vertex_z = config_.vertex_z_m;

    for (int i = 0; i < num_points; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(num_points - 1);
        const double r = alpha * aperture;
        const double z = vertex_z + (r * r) / (4.0 * focal_length);
        points.push_back({r, z});
    }

    return points;
}
