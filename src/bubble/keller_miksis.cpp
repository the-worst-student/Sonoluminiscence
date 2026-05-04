#include "bubble/keller_miksis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double Cube(double value) {
    return value * value * value;
}

}  // namespace

KellerMiksisModel::KellerMiksisModel(const BubblePhysicalParameters& parameters, const BubbleExcitationInput& excitation)
    : parameters_(parameters),
      excitation_(excitation),
      equilibrium_radius_m_(excitation.EquilibriumRadiusM),
      hard_core_radius_m_(parameters.HardCoreRatio * excitation.EquilibriumRadiusM) {
    if (equilibrium_radius_m_ <= 0.0) {
        throw std::runtime_error("Equilibrium bubble radius must be positive.");
    }

    if (hard_core_radius_m_ <= 0.0) {
        throw std::runtime_error("Hard-core radius must be positive.");
    }

    if (hard_core_radius_m_ >= equilibrium_radius_m_) {
        throw std::runtime_error("Hard-core radius must be smaller than equilibrium radius.");
    }

    if (excitation.LiquidTemperatureK <= 0.0) {
        throw std::runtime_error("Liquid temperature must be positive.");
    }

    const double initial_gas_pressure_pa = excitation.StaticPressurePa + 2.0 * parameters_.LiquidSurfaceTensionNPerM / equilibrium_radius_m_;
    const double initial_effective_volume_m3 = EffectiveVolumeM3(equilibrium_radius_m_);

    argon_amount_mol_ = initial_gas_pressure_pa * initial_effective_volume_m3 /
                        (parameters_.GasMolarConstantJPerMolK * excitation.LiquidTemperatureK);

    if (argon_amount_mol_ <= 0.0) {
        throw std::runtime_error("Computed argon amount must be positive.");
    }
}

KmBubbleStateDerivative KellerMiksisModel::Evaluate(double t, const KmBubbleState& state) const {
    if (state.R <= 1.001 * hard_core_radius_m_) {
        throw std::runtime_error("Bubble radius is too close to hard-core radius.");
    }

    if (state.Tg <= 0.0) {
        throw std::runtime_error("Gas temperature must be positive.");
    }

    const double R = state.R;
    const double U = state.U;
    const double Tg = state.Tg;

    const double c = parameters_.LiquidSoundSpeedMPerS;
    const double mach = U / c;

    if (std::abs(U) > 0.8 * c) {
        throw std::runtime_error("Bubble wall velocity exceeded liquid Mach safety limit.");
    }

    if (std::abs(1.0 - mach) < 0.05) {
        throw std::runtime_error("Keller-Miksis compressibility denominator is too close to zero.");
    }

    const double pg = GasPressurePa(state);

    if (!std::isfinite(pg) || pg <= 0.0) {
        throw std::runtime_error("Gas pressure must be positive and finite.");
    }

    const double pinf = ExternalPressurePa(t);
    const double dpinf_dt = ExternalPressureDerivativePaS(t);

    const double v_eff = EffectiveVolumeM3(R);
    const double dv_eff_dt = EffectiveVolumeDerivativeM3S(state);
    const double q_out = HeatLossRateW(state);

    const double dTg_dt = (-pg * dv_eff_dt - q_out) /
                          (argon_amount_mol_ * parameters_.GasIsochoricHeatCapacityJPerMolK);

    const double dpg_dt = pg * (dTg_dt / Tg - dv_eff_dt / v_eff);

    const double delta_p = pg - pinf -
                           2.0 * parameters_.LiquidSurfaceTensionNPerM / R -
                           4.0 * parameters_.LiquidDynamicViscosityPaS * U / R;

    const double d_pg_minus_pinf_dt = dpg_dt - dpinf_dt;

    const double rho = parameters_.LiquidDensityKgPerM3;

    const double denominator = R * (1.0 - mach);

    if (std::abs(denominator) < 1.0e-30) {
        throw std::runtime_error("Keller-Miksis denominator is too close to zero.");
    }

    const double numerator = (1.0 / rho) * (1.0 + mach) * delta_p +
                             R / (rho * c) * d_pg_minus_pinf_dt -
                             1.5 * U * U * (1.0 - U / (3.0 * c));

    KmBubbleStateDerivative derivative;
    derivative.dRdt = U;
    derivative.dUdt = numerator / denominator;
    derivative.dTgdt = dTg_dt;

    return derivative;
}

double KellerMiksisModel::GasPressurePa(const KmBubbleState& state) const {
    if (state.Tg <= 0.0) {
        throw std::runtime_error("Gas temperature must be positive.");
    }

    const double effective_volume_m3 = EffectiveVolumeM3(state.R);
    const double pressure_pa = argon_amount_mol_ * parameters_.GasMolarConstantJPerMolK * state.Tg / effective_volume_m3;

    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
        throw std::runtime_error("Gas pressure must be positive and finite.");
    }

    return pressure_pa;
}

double KellerMiksisModel::ExternalPressurePa(double t) const {
    return excitation_.StaticPressurePa +
           excitation_.DriveAmplitudePa * std::cos(excitation_.OmegaRadS * t + excitation_.DrivePhaseRad);
}

double KellerMiksisModel::ExternalPressureDerivativePaS(double t) const {
    return -excitation_.DriveAmplitudePa * excitation_.OmegaRadS *
           std::sin(excitation_.OmegaRadS * t + excitation_.DrivePhaseRad);
}

double KellerMiksisModel::ThermalLayerThicknessM(const KmBubbleState& state) const {
    constexpr double kVelocityEpsilonMPerS = 1.0e-9;

    if (state.R <= 0.0) {
        throw std::runtime_error("Bubble radius must be positive.");
    }

    if (excitation_.OmegaRadS <= 0.0) {
        throw std::runtime_error("Angular frequency must be positive.");
    }

    const double acoustic_time_s = 1.0 / excitation_.OmegaRadS;
    const double collapse_time_s = state.R / (std::abs(state.U) + kVelocityEpsilonMPerS);
    const double local_time_s = std::min(acoustic_time_s, collapse_time_s);

    const double diffusive_layer_m =
        parameters_.ThermalLayerCoefficient * std::sqrt(parameters_.GasThermalDiffusivityM2PerS * local_time_s);

    return std::min(state.R, diffusive_layer_m);
}

double KellerMiksisModel::HeatLossRateW(const KmBubbleState& state) const {
    const double delta_t_m = ThermalLayerThicknessM(state);

    if (delta_t_m <= 0.0) {
        throw std::runtime_error("Thermal layer thickness must be positive.");
    }

    return 4.0 * kPi * state.R * state.R * parameters_.GasThermalConductivityWPerMK *
           (state.Tg - excitation_.LiquidTemperatureK) / delta_t_m;
}

double KellerMiksisModel::HardCoreRadiusM() const {
    return hard_core_radius_m_;
}

double KellerMiksisModel::ArgonAmountMol() const {
    return argon_amount_mol_;
}

double KellerMiksisModel::EffectiveVolumeM3(double radius_m) const {
    if (radius_m <= hard_core_radius_m_) {
        throw std::runtime_error("Radius must be larger than hard-core radius.");
    }

    const double volume_m3 = 4.0 * kPi / 3.0 * (Cube(radius_m) - Cube(hard_core_radius_m_));

    if (!std::isfinite(volume_m3) || volume_m3 <= 0.0) {
        throw std::runtime_error("Effective gas volume must be positive and finite.");
    }

    return volume_m3;
}

double KellerMiksisModel::EffectiveVolumeDerivativeM3S(const KmBubbleState& state) const {
    return 4.0 * kPi * state.R * state.R * state.U;
}