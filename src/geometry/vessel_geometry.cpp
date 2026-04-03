#include "geometry/vessel_geometry.hpp"

VesselGeometry::VesselGeometry(const VesselConfig& config) : config_(config) {}

VesselGeometryData VesselGeometry::BuildData() const {
    VesselGeometryData data;
    data.radius_m = config_.radius_m;
    data.height_m = config_.height_m;
    return data;
}

