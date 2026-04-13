#pragma once

class GasModel {
public:
    static double EffectiveVolume(
        double radius_m,
        double hard_core_radius_m);

    static double PressurePa(
        double gas_moles_mol,
        double gas_temperature_k,
        double radius_m,
        double hard_core_radius_m);

    static double InitialMolesFromEquilibrium(
        double equilibrium_radius_m,
        double hard_core_radius_m,
        double gas_temperature_k,
        double static_pressure_pa,
        double surface_tension_n_m);
};