#include <bubble/rayleigh_plesset.hpp>
#include <stdexcept>

double RayleighPlesset::RadialAccelerationMPerS2(
    double radius_m,
    double radius_velocity_m_s,
    double liquid_density_kg_m3,
    double liquid_viscosity_pa_s,
    double surface_tension_n_m,
    double gas_pressure_pa,
    double external_pressure_pa) {
    if (radius_m <= 0.0) {
        throw std::runtime_error("Bubble radius must be positive");
    }
    if (liquid_density_kg_m3 <= 0.0) {
        throw std::runtime_error("Liquid density must be positive");
    }
    const double pressure_term =
        (gas_pressure_pa - external_pressure_pa -
         2.0 * surface_tension_n_m / radius_m -
         4.0 * liquid_viscosity_pa_s * radius_velocity_m_s / radius_m) /
        liquid_density_kg_m3;
    return (pressure_term - 1.5 * radius_velocity_m_s * radius_velocity_m_s) /
           radius_m;
}

