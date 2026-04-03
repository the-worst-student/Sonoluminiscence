#pragma once

#include <utility>
#include <vector>

#include <core/config.hpp>

struct ReflectorPoint {
    double r_m;
    double z_m;
};

class ReflectorGeometry {
public:
    explicit ReflectorGeometry(const ReflectorConfig& config);

    std::vector<ReflectorPoint> BuildProfile(int num_points) const;
private:
    ReflectorConfig config_;
};

