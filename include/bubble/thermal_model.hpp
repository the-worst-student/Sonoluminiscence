#pragma once

class ThermalModel {
public:
    static double HeatLossRateW(
        double radius_m,
        double gas_temperature_k,
        double liquid_temperature,
        double thermal_conductivity_w_m_k,
        double thermal_layer_thickness_m);
};