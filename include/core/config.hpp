#pragma once

#include <string>
#include <vector>

struct ProjectInfoConfig {
    std::string name;
    std::string run_id;
};

struct VesselConfig {
    std::string type;
    double radius_m;
    double height_m;
};

struct ReflectorConfig {
    bool enabled;
    std::string type;
    double focal_length_m;
    double vertex_z_m;
    double aperture_radius_m;
};

struct SourcePatchConfig {
    std::string type;
    double r_min_m;
    double r_max_m;
    double z_m;
};

struct BubblePositionConfig {
    double r_m;
    double z_m;
};

struct GeometryConfig {
    std::string formulation;
    VesselConfig vessel;
    ReflectorConfig reflector;
    SourcePatchConfig source_patch;
    BubblePositionConfig bubble_position;
};

struct MeshConfig {
    double h_bulk_m;
    double h_source_m;
    double h_reflector_m;
    double h_bubble_zone_m;
};

struct LiquidConfig {
    double density_kg_m3;
    double sound_speed_m_s;
    double viscosity_pa_s;
    double surface_tension_n_m;
    double static_pressure_pa;
    double temperature_k;
};

struct AcousticSourceConfig {
    int boundary_tag;
    double normal_velocity_m_s;
    std::string mode;
};

struct AcousticsOutputConfig {
    bool sample_bubble_point;
    bool compute_focus_region;
};

struct AcousticsConfig {
    std::string formulation;
    double frequency_hz;
    std::string time_convention;
    AcousticSourceConfig source;
    std::vector<int> rigid_boundaries;
    AcousticsOutputConfig output;
};

struct BubbleDriveConfig {
    double pressure_amplitude_pa;
    double phase_rad;
};

struct BubbleGasConfig {
    double hard_core_radius_m;
    double heat_capacity_cv_j_mol_k;
    double thermal_conductivity_w_m_k;
    double thermal_layer_thickness_m;
};

struct BubbleIntegrationConfig {
    double time_step_s;
    double final_time_s;
    int output_every;
};

struct BubbleConfig {
    double equilibrium_radius_m;
    double initial_radius_m;
    double initial_velocity_m_s;
    double initial_temperature_k;
    BubbleDriveConfig drive;
    BubbleGasConfig gas;
    BubbleIntegrationConfig integration;
};

struct ProjectConfig {
    ProjectInfoConfig project;
    GeometryConfig geometry;
    MeshConfig mesh{};
    LiquidConfig liquid{};
    AcousticsConfig acoustics;
    BubbleConfig bubble{};
};
