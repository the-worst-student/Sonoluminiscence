#include <bubble/bubble_rhs.hpp>
#include <bubble/gas_model.hpp>
#include <bubble/rayleigh_plesset.hpp>
#include <bubble/thermal_model.hpp>

namespace {
constexpr double cPi = 3.14159265358979323846;
}

BubbleRhs::BubbleRhs(const ProjectConfig& config) : Parameters((BuildParameters(config))),
    ExternalPressure(
        Parameters.StaticPressurePa,
        config.bubble.drive.pressure_amplitude_pa,
        Parameters.AngularFrequencyRadS,
        config.bubble.drive.phase_rad) {};

BubbleDerivative BubbleRhs::Evaluate(const double time_s, const BubbleState& state) const {
    const double external_pressure = ComputeExternalPressure(time_s);
    const double gas_pressure = ComputeGasPressure(state);

    const double radius_rate = state.RadiusVelocityMPerS;
    const double radius_acceleration = RayleighPlesset::RadialAccelerationMPerS2(
        state.RadiusM,
        state.RadiusVelocityMPerS,
        Parameters.LiquidDensityKgPerM3,
        Parameters.LiquidViscosityPaS,
        Parameters.SurfaceTensionNPerM,
        gas_pressure,
        external_pressure);

    const double volume_rate = 4.0 * cPi * state.RadiusM * state.RadiusM * state.RadiusVelocityMPerS;
    const double heat_loss =
    ThermalModel::HeatLossRateW(
        state.RadiusM,
        state.GasTemperatureK,
        Parameters.LiquidTemperatureK,
        Parameters.GasThermalConductivityWPerMK,
        Parameters.ThermalLayerThicknessM);

    const double gas_temperature_rate = (-gas_pressure * volume_rate - heat_loss) / (Parameters.GasMolesMol * Parameters.GasCvJPerMolK);

    BubbleDerivative derivative;
    derivative.RadiusRateMPerS = radius_rate;
    derivative.RadiusAccelerationMPerS2 = radius_acceleration;
    derivative.GasTemperatureRateKPerS = gas_temperature_rate;

    return derivative;
}

double BubbleRhs::ComputeExternalPressure(const double time_s) const {
    return ExternalPressure.Evaluate(time_s);
}

double BubbleRhs::ComputeGasPressure(const BubbleState& state) const {
    return GasModel::PressurePa(
        Parameters.GasMolesMol,
        state.GasTemperatureK,
        state.RadiusM,
        Parameters.HardCoreRadiusM);
}

BubbleState BubbleRhs::BuildInitialState() const {
    BubbleState state;
    state.RadiusM = Parameters.InitialRadiusM;
    state.RadiusVelocityMPerS = Parameters.InitialVelocityMPerS;
    state.GasTemperatureK = Parameters.InitialTemperatureK;
    return state;
}

const BubbleModelParameters& BubbleRhs::GetParameters() const {
    return Parameters;
}

BubbleModelParameters BubbleRhs::BuildParameters(const ProjectConfig& config) {
    BubbleModelParameters parameters;
    parameters.LiquidDensityKgPerM3 = config.liquid.density_kg_m3;
    parameters.LiquidViscosityPaS = config.liquid.viscosity_pa_s;
    parameters.SurfaceTensionNPerM = config.liquid.surface_tension_n_m;
    parameters.StaticPressurePa = config.liquid.static_pressure_pa;
    parameters.LiquidTemperatureK = config.liquid.temperature_k;
    parameters.AngularFrequencyRadS = 2.0 * cPi * config.acoustics.frequency_hz;
    parameters.EquilibriumRadiusM = config.bubble.equilibrium_radius_m;
    parameters.InitialRadiusM = config.bubble.initial_radius_m;
    parameters.InitialVelocityMPerS = config.bubble.initial_velocity_m_s;
    parameters.InitialTemperatureK = config.bubble.initial_temperature_k;
    parameters.HardCoreRadiusM = config.bubble.gas.hard_core_radius_m;
    parameters.GasCvJPerMolK = config.bubble.gas.heat_capacity_cv_j_mol_k;
    parameters.GasThermalConductivityWPerMK = config.bubble.gas.thermal_conductivity_w_m_k;
    parameters.ThermalLayerThicknessM = config.bubble.gas.thermal_layer_thickness_m;
    parameters.GasMolesMol =
    GasModel::InitialMolesFromEquilibrium(
        parameters.EquilibriumRadiusM,
        parameters.HardCoreRadiusM,
        parameters.InitialTemperatureK,
        parameters.StaticPressurePa,
        parameters.SurfaceTensionNPerM);

    return parameters;
}