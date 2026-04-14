#pragma once

class RayleighPlesset {
public:
    static double RadialAccelerationMPerS2(
        double radius_m,
        double radius_velocity_m_s,
        double liquid_density_kg_m3,
        double liquid_viscosity_pa_s,
        double surface_tension_n_m,
        double gas_pressure_pa,
        double external_pressure_pa);
};