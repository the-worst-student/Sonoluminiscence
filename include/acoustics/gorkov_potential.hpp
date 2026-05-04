#pragma once

#include <acoustics/field_sampler.hpp>
#include <core/config.hpp>

struct GorkovParameters {
    double BubbleRadiusM = 0.0;
    double LiquidDensityKgM3 = 0.0;
    double SoundSpeedMS = 0.0;
    double AngularFrequencyRadS = 0.0;
    double CompressibilityContrast = 1.0;
    double DensityContrast = -2.0;
};

class GorkovPotential {
public:
    static GorkovParameters FromConfig(
        const ProjectConfig& config,
        double angular_frequency_rad_s);

    static void AddGorkovFields(
        AcousticFieldData* field_data,
        const GorkovParameters& parameters);
};
