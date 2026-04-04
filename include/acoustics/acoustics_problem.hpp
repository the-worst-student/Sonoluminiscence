#pragma once

#include <complex>
#include <vector>

#include <core/config.hpp>

enum class TimeConvention {
  ExpMinusIwt,
  ExpPlusIwt,
};

struct AcousticsProblem {
    double FrequencyHz = 0.0;
    double Omega = 0.0;
    double Density = 0.0;
    double SoundSpeed = 0.0;
    double WaveNumber = 0.0;
    double SourceNormalVelocity;

    TimeConvention Convention = TimeConvention::ExpMinusIwt;

    int SourceBoundaryTag = -1;
    std::vector<int> RigidBoundaryTags;
    std::complex<double> SourceFlux = {0.0, 0.0};

    static AcousticsProblem FromConfig(const ProjectConfig& config);
    void Validate() const;

    bool IsRigidBoundary(int tag) const;
    bool UseMinusIwmConvention() const;
};