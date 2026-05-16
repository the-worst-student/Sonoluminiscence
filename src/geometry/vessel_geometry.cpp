#include "geometry/vessel_geometry.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace {

constexpr double kPi = 3.14159265358979323846;

double LinearRadius(double bottom_radius_m, double top_radius_m, double alpha) {
    return bottom_radius_m + alpha * (top_radius_m - bottom_radius_m);
}

}  // namespace

VesselGeometry::VesselGeometry(const VesselConfig& config) : config_(config) {}

VesselGeometryData VesselGeometry::BuildData() const {
    ValidateConfig();

    VesselGeometryData data;
    data.radius_m = config_.radius_m;
    data.height_m = config_.height_m;

    const int point_count = ProfilePointCount();

    data.outer_wall.reserve(static_cast<std::size_t>(point_count));

    for (int i = 0; i < point_count; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(point_count - 1);
        const double z_m = alpha * config_.height_m;
        const double r_m = RadiusAtAlpha(alpha);

        if (r_m <= 0.0) {
            throw std::invalid_argument("Vessel wall radius must be positive");
        }

        data.outer_wall.push_back({r_m, z_m});
    }

    return data;
}

double VesselGeometry::RadiusAtZ(double z_m) const {
    ValidateConfig();

    if (z_m < 0.0 || z_m > config_.height_m) {
        throw std::invalid_argument("z-coordinate is outside vessel height");
    }

    const double alpha = z_m / config_.height_m;
    return RadiusAtAlpha(alpha);
}

void VesselGeometry::ValidateConfig() const {
    if (config_.radius_m <= 0.0) {
        throw std::invalid_argument("Vessel radius_m must be positive");
    }

    if (config_.height_m <= 0.0) {
        throw std::invalid_argument("Vessel height_m must be positive");
    }

    if (config_.bottom_radius_m <= 0.0) {
        throw std::invalid_argument("Vessel bottom_radius_m must be positive");
    }

    if (config_.top_radius_m <= 0.0) {
        throw std::invalid_argument("Vessel top_radius_m must be positive");
    }

    if (config_.profile_points < 2) {
        throw std::invalid_argument("Vessel profile_points must be at least 2");
    }

    if (config_.bulge_m < 0.0) {
        throw std::invalid_argument("Vessel bulge_m must be non-negative");
    }

    if (config_.neck_m < 0.0) {
        throw std::invalid_argument("Vessel neck_m must be non-negative");
    }

    if (config_.type != "cylinder" &&
        config_.type != "conical" &&
        config_.type != "barrel" &&
        config_.type != "hourglass") {
        throw std::invalid_argument("Unsupported vessel type: " + config_.type);
    }
}

int VesselGeometry::ProfilePointCount() const {
    if (config_.type == "cylinder" || config_.type == "conical") {
        return 2;
    }

    return config_.profile_points;
}

double VesselGeometry::RadiusAtAlpha(double alpha) const {
    const double linear_radius =
        LinearRadius(config_.bottom_radius_m, config_.top_radius_m, alpha);

    if (config_.type == "cylinder") {
        return config_.radius_m;
    }

    if (config_.type == "conical") {
        return linear_radius;
    }

    if (config_.type == "barrel") {
        return linear_radius + config_.bulge_m * std::sin(kPi * alpha);
    }

    if (config_.type == "hourglass") {
        return linear_radius - config_.neck_m * std::sin(kPi * alpha);
    }

    throw std::invalid_argument("Unsupported vessel type: " + config_.type);
}