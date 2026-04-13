#include "bubble/thermal_model.hpp"

#include <stdexcept>

namespace {
constexpr double cPi   = 3.14159265358979323846;
}

double ThermalModel::HeatLossRateW(
    const double radius_m,
    const double gas_temperature_k,
    const double liquid_temperature_k,
    const double thermal_conductivity_w_m_k,
    const double thermal_layer_thickness_m) {
    if (radius_m <= 0.0) {
        throw std::runtime_error("Bubble radius must be positive");
    }

    if (thermal_layer_thickness_m <= 0.0) {
        throw std::runtime_error("Thermal layer thickness must be positive");
    }

    return 4.0 * cPi   * radius_m * radius_m *
           thermal_conductivity_w_m_k *
           (gas_temperature_k - liquid_temperature_k) /
           thermal_layer_thickness_m;
}