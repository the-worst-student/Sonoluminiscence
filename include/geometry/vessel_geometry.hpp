#pragma once

#include "core/config.hpp"

struct VesselGeometryData {
    double radius_m;
    double height_m;
};

class VesselGeometry {
public:
    explicit VesselGeometry(const VesselConfig& config);

    VesselGeometryData BuildData() const;

private:
    VesselConfig config_;
};