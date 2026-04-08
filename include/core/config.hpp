#pragma once
#include <vector>
#include <string>

struct ProjectInfoConfig {
    std::string name;
    std::string run_id;
};

struct VesselConfig {
    std::string type;
    double radius_m;
    double height_m;
};

// Reflector's form equation = z0 + r^2/4f
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

struct ProjectConfig {
    ProjectInfoConfig project;
    GeometryConfig geometry;
    MeshConfig mesh{};
    LiquidConfig liquid{};
    AcousticsConfig acoustics;
};