#include "coupling/pressure_coupling.hpp"

#include <cmath>

PressureCoupling::PressureCoupling(
    const double static_pressure_pa,
    const double acoustic_amplitude_pa,
    const double angular_frequency_rad_s,
    const double phase_rad)
    : StaticPressurePa(static_pressure_pa),
      AcousticAmplitudePa(acoustic_amplitude_pa),
      AngularFrequencyRadS(angular_frequency_rad_s),
      PhaseRad(phase_rad) {}

double PressureCoupling::AcousticComponent(const double time_s) const {
    return AcousticAmplitudePa * std::cos(AngularFrequencyRadS * time_s + PhaseRad);
}

double PressureCoupling::Evaluate(const double time_s) const {
    return StaticPressurePa + AcousticComponent(time_s);
}