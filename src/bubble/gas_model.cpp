#include <bubble/gas_model.hpp>
#include <cmath>
#include <stdexcept>

namespace {
constexpr double cPi   = 3.14159265358979323846;
constexpr double cUniversalGasConstant                  = 8.31446261815324;
}

double GasModel::EffectiveVolume(const double radius_m, const double hard_core_radius_m) {
    if (radius_m <= hard_core_radius_m) {
        throw std::runtime_error("Bubble radius must be greater than hard-core radius");
    }
    return (4.0 * cPi / 3.0) * (radius_m * radius_m * radius_m - hard_core_radius_m * hard_core_radius_m * hard_core_radius_m);
}

double GasModel::PressurePa(double gas_moles_mol, double gas_temperature_k, double radius_m, double hard_core_radius_m) {
    if (gas_moles_mol <= 0.0) {
        throw std::runtime_error("Gas amount must be positive");
    }
    const double effective_volume = EffectiveVolume(radius_m, hard_core_radius_m);
    return gas_moles_mol * cUniversalGasConstant                  * gas_temperature_k / effective_volume;
}

double GasModel::InitialMolesFromEquilibrium(double equilibrium_radius_m, double hard_core_radius_m,
    double gas_temperature_k, double static_pressure_pa, double surface_tension_n_m) {
    if (equilibrium_radius_m <= 0.0) {
        throw std::runtime_error("Equilibrium radius must be positive");
    }
    const double gas_pressure_equilibrium = static_pressure_pa + 2.0 * surface_tension_n_m / equilibrium_radius_m;
    const double effective_volume = EffectiveVolume(equilibrium_radius_m, hard_core_radius_m);
    return gas_pressure_equilibrium * effective_volume / (cUniversalGasConstant                   * gas_temperature_k);

}