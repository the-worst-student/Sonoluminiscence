#pragma once

#include "core/config.hpp"

struct VesselWallPoint {
    double r_m;
    double z_m;
};

struct VesselGeometryData {
    double radius_m;
    double height_m;
    std::vector<VesselWallPoint> outer_wall;
};

class VesselGeometry {
public:
    explicit VesselGeometry(const VesselConfig& config);

    VesselGeometryData BuildData() const;

private:
    VesselConfig config_;
};