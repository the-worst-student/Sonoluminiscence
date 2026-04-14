#pragma once

#include <complex>

class PressureCoupling {
public:
    PressureCoupling(
        double static_pressure_pa,
        double acoustic_amplitude_pa,
        double angular_frequency_rad_s,
        double phase_rad);

    static PressureCoupling FromComplexPressure(
        double static_pressure_pa,
        double angular_frequency_rad_s,
        std::complex<double> pb);

    double AcousticComponent(double time_s) const;
    double Evaluate(double time_s) const;

    double GetStaticPressurePa() const;
    double GetAcousticAmplitudePa() const;
    double GetAngularFrequencyRadS() const;
    double GetPhaseRad() const;

private:
    double StaticPressurePa = 0.0;
    double AcousticAmplitudePa = 0.0;
    double AngularFrequencyRadS = 0.0;
    double PhaseRad = 0.0;
};