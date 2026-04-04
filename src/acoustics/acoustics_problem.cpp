#include "acoustics/acoustics_problem.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

TimeConvention ParseTimeConvention(const std::string& value) {
  if (value == "exp_minus_iwt") {
    return TimeConvention::ExpMinusIwt;
  }
  if (value == "exp_plus_iwt") {
    return TimeConvention::ExpPlusIwt;
  }

  throw std::invalid_argument("Unsupported acoustics.time_convention: " + value);
}

}  // namespace

AcousticsProblem AcousticsProblem::FromConfig(const ProjectConfig& config) {
  AcousticsProblem problem;

  problem.FrequencyHz = config.acoustics.frequency_hz;
  problem.Density = config.liquid.density_kg_m3;
  problem.SoundSpeed = config.liquid.sound_speed_m_s;
  problem.SourceNormalVelocity = config.acoustics.source.normal_velocity_m_s;
  problem.Convention = ParseTimeConvention(config.acoustics.time_convention);
  problem.SourceBoundaryTag = config.acoustics.source.boundary_tag;
  problem.RigidBoundaryTags = config.acoustics.rigid_boundaries;

  constexpr double kPi = 3.14159265358979323846;
  problem.Omega = 2.0 * kPi * problem.FrequencyHz;
  problem.WaveNumber = problem.Omega / problem.SoundSpeed;

  const double source_flux_imag =
      problem.Omega * problem.Density * problem.SourceNormalVelocity;

  if (problem.Convention == TimeConvention::ExpMinusIwt) {
    problem.SourceFlux = std::complex<double>(0.0, source_flux_imag);
  } else {
    problem.SourceFlux = std::complex<double>(0.0, -source_flux_imag);
  }

  problem.Validate();
  return problem;
}

void AcousticsProblem::Validate() const {
  if (FrequencyHz <= 0.0) {
    throw std::invalid_argument("Acoustics frequency must be positive");
  }
  if (Density <= 0.0) {
    throw std::invalid_argument("Liquid density must be positive");
  }
  if (SoundSpeed <= 0.0) {
    throw std::invalid_argument("Liquid sound speed must be positive");
  }
  if (WaveNumber <= 0.0) {
    throw std::invalid_argument("Wave number must be positive");
  }
  if (SourceBoundaryTag <= 0) {
    throw std::invalid_argument("Source boundary tag must be positive");
  }
  if (std::find(
          RigidBoundaryTags.begin(),
          RigidBoundaryTags.end(),
          SourceBoundaryTag) != RigidBoundaryTags.end()) {
    throw std::invalid_argument(
        "Source boundary tag must not belong to rigid boundaries");
  }
}

bool AcousticsProblem::IsRigidBoundary(const int tag) const {
  return std::find(RigidBoundaryTags.begin(), RigidBoundaryTags.end(), tag) !=
         RigidBoundaryTags.end();
}

bool AcousticsProblem::UsesMinusIwtConvention() const {
  return Convention == TimeConvention::ExpMinusIwt;
}
