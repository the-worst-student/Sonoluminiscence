#include "bubble/adaptive_ode_solver.hpp"
#include "bubble/bubble_simulation_result.hpp"
#include "bubble/keller_miksis.hpp"
#include "io/bubble_excitation_reader.hpp"
#include "io/yaml_reader.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kSimulationCycles = 10.0;

struct CandidateSimulationOutput {
    BubbleSimulationResult Result;
    std::vector<BubbleSample> Samples;
};

BubblePhysicalParameters MakePhysicalParametersFromConfig(const ProjectConfig& config) {
    BubblePhysicalParameters parameters;

    parameters.LiquidDensityKgPerM3 = config.liquid.density_kg_m3;
    parameters.LiquidSoundSpeedMPerS = config.liquid.sound_speed_m_s;
    parameters.LiquidDynamicViscosityPaS = config.liquid.viscosity_pa_s;
    parameters.LiquidSurfaceTensionNPerM = config.liquid.surface_tension_n_m;

    parameters.GasIsochoricHeatCapacityJPerMolK = config.bubble.gas.heat_capacity_cv_j_mol_k;
    parameters.GasThermalConductivityWPerMK = config.bubble.gas.thermal_conductivity_w_m_k;

    if (config.bubble.equilibrium_radius_m > 0.0 && config.bubble.gas.hard_core_radius_m > 0.0) {
        parameters.HardCoreRatio = config.bubble.gas.hard_core_radius_m / config.bubble.equilibrium_radius_m;
    }

    return parameters;
}

BubbleState MakeInitialBubbleState(const BubbleExcitationInput& excitation) {
    BubbleState state;
    state.R = excitation.EquilibriumRadiusM;
    state.U = 0.0;
    state.Tg = excitation.LiquidTemperatureK;
    return state;
}

AdaptiveOdeOptions MakeDefaultOdeOptions(double drive_period_s) {
    AdaptiveOdeOptions options;
    options.InitialStepS = drive_period_s / 10000.0;
    options.MaxStepS = drive_period_s / 1000.0;
    options.MinStepS = drive_period_s / 1.0e9;
    options.RelativeTolerance = 1.0e-6;
    options.AbsoluteToleranceR = 1.0e-12;
    options.AbsoluteToleranceU = 1.0e-3;
    options.AbsoluteToleranceT = 1.0e-2;
    options.MaxSteps = 1000000;
    return options;
}

BubbleSimulationResult MakeNumericallyInvalidResult(const BubbleExcitationInput& excitation, const BubblePhysicalParameters& parameters, const std::string& reason) {
    BubbleSimulationResult result;

    result.CaseId = excitation.CaseId;
    result.Status = "NumericallyInvalid";
    result.FailureReason = reason;

    result.FrequencyHz = excitation.FrequencyHz;
    result.OmegaRadS = excitation.OmegaRadS;
    result.DriveAmplitudePa = excitation.DriveAmplitudePa;
    result.DrivePhaseRad = excitation.DrivePhaseRad;

    result.StaticPressurePa = excitation.StaticPressurePa;
    result.LiquidTemperatureK = excitation.LiquidTemperatureK;

    result.R0M = excitation.EquilibriumRadiusM;
    result.HardCoreRadiusM = parameters.HardCoreRatio * excitation.EquilibriumRadiusM;

    return result;
}

CandidateSimulationOutput MakeNumericallyInvalidOutput(const BubbleExcitationInput& excitation, const BubblePhysicalParameters& parameters, const std::string& reason) {
    CandidateSimulationOutput output;
    output.Result = MakeNumericallyInvalidResult(excitation, parameters, reason);
    return output;
}

CandidateSimulationOutput SimulateBubbleCandidate(const BubbleExcitationInput& excitation, const BubblePhysicalParameters& parameters) {
    try {
        if (excitation.OmegaRadS <= 0.0) {
            return MakeNumericallyInvalidOutput(excitation, parameters, "NonPositiveAngularFrequency");
        }

        KellerMiksisModel model(parameters, excitation);

        const BubbleState initial_state = MakeInitialBubbleState(excitation);
        const double drive_period_s = 2.0 * kPi / excitation.OmegaRadS;

        AdaptiveOdeOptions ode_options = MakeDefaultOdeOptions(drive_period_s);
        AdaptiveOdeSolver solver;

        const auto rhs = [&model](double t, const BubbleState& state) {
            return model.Evaluate(t, state);
        };

        const double t_start = 0.0;
        const double t_end = kSimulationCycles * drive_period_s;

        const AdaptiveOdeResult ode_result = solver.Solve(initial_state, t_start, t_end, ode_options, rhs);

        CandidateSimulationOutput output;
        output.Samples = ode_result.Samples;
        output.Result = AnalyzeBubbleSimulationResult(excitation, model, parameters, ode_result.Samples, ode_result.Success, ode_result.FailureReason);
        return output;
    } catch (const std::exception& exception) {
        return MakeNumericallyInvalidOutput(excitation, parameters, exception.what());
    }
}

int BoolToInt(bool value) {
    return value ? 1 : 0;
}

std::string EscapeCsvField(const std::string& value) {
    bool needs_quotes = false;

    for (char symbol : value) {
        if (symbol == ',' || symbol == '"' || symbol == '\n' || symbol == '\r') {
            needs_quotes = true;
            break;
        }
    }

    if (!needs_quotes) {
        return value;
    }

    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');

    for (char symbol : value) {
        if (symbol == '"') {
            escaped += "\"\"";
        } else {
            escaped.push_back(symbol);
        }
    }

    escaped.push_back('"');
    return escaped;
}

void WriteBubbleResultsCsv(const std::string& output_csv_path, const std::vector<BubbleSimulationResult>& results) {
    const std::filesystem::path path(output_csv_path);

    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream file(output_csv_path);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output CSV: " + output_csv_path);
    }

    file << "case_id,status,failure_reason,"
         << "frequency_hz,omega_rad_s,drive_amplitude_pa,drive_phase_rad,"
         << "static_pressure_pa,liquid_temperature_k,"
         << "R0_m,h_m,R_min_m,R_max_m,"
         << "T_max_k,p_g_max_pa,U_max_m_s,"
         << "K_R,K_exp,K_range,M_l,K_h,"
         << "min_radius_detected,mechanical_collapse_ok,thermal_heating_ok,"
         << "thermomechanical_luminescence_candidate\n";

    for (const BubbleSimulationResult& result : results) {
        file << EscapeCsvField(result.CaseId) << ','
             << EscapeCsvField(result.Status) << ','
             << EscapeCsvField(result.FailureReason) << ','
             << result.FrequencyHz << ','
             << result.OmegaRadS << ','
             << result.DriveAmplitudePa << ','
             << result.DrivePhaseRad << ','
             << result.StaticPressurePa << ','
             << result.LiquidTemperatureK << ','
             << result.R0M << ','
             << result.HardCoreRadiusM << ','
             << result.RMinM << ','
             << result.RMaxM << ','
             << result.TMaxK << ','
             << result.PgMaxPa << ','
             << result.UMaxMS << ','
             << result.CompressionRatio << ','
             << result.ExpansionRatio << ','
             << result.DynamicRange << ','
             << result.LiquidMach << ','
             << result.HardCoreRatio << ','
             << BoolToInt(result.MinRadiusDetected) << ','
             << BoolToInt(result.MechanicalCollapseOk) << ','
             << BoolToInt(result.ThermalHeatingOk) << ','
             << BoolToInt(result.ThermoMechanicalLuminescenceCandidate)
             << '\n';
    }
}

bool IsBetterSbsCandidate(const BubbleSimulationResult& lhs, const BubbleSimulationResult& rhs) {
    if (lhs.ThermoMechanicalLuminescenceCandidate != rhs.ThermoMechanicalLuminescenceCandidate) {
        return lhs.ThermoMechanicalLuminescenceCandidate > rhs.ThermoMechanicalLuminescenceCandidate;
    }

    if (lhs.TMaxK != rhs.TMaxK) {
        return lhs.TMaxK > rhs.TMaxK;
    }

    if (lhs.CompressionRatio != rhs.CompressionRatio) {
        return lhs.CompressionRatio > rhs.CompressionRatio;
    }

    if (lhs.PgMaxPa != rhs.PgMaxPa) {
        return lhs.PgMaxPa > rhs.PgMaxPa;
    }

    if (lhs.LiquidMach != rhs.LiquidMach) {
        return lhs.LiquidMach > rhs.LiquidMach;
    }

    return lhs.CaseId < rhs.CaseId;
}

std::string MakeBestCandidatesPath(const std::string& output_csv_path) {
    const std::filesystem::path output_path(output_csv_path);
    const std::filesystem::path parent_path = output_path.parent_path();

    if (parent_path.empty()) {
        return "best_sbs_candidates.csv";
    }

    return (parent_path / "best_sbs_candidates.csv").string();
}

void WriteBestSbsCandidatesCsv(const std::string& best_csv_path,
                               const std::vector<BubbleSimulationResult>& results,
                               std::size_t max_best_count) {
    std::vector<BubbleSimulationResult> sorted_results = results;

    std::sort(sorted_results.begin(), sorted_results.end(), IsBetterSbsCandidate);

    if (sorted_results.size() > max_best_count) {
        sorted_results.resize(max_best_count);
    }

    WriteBubbleResultsCsv(best_csv_path, sorted_results);
}

void WriteBubbleTimeSeriesCsv(const std::string& timeseries_dir, const BubbleExcitationInput& excitation, const BubblePhysicalParameters& parameters, const std::vector<BubbleSample>& samples) {
    if (samples.empty()) {
        return;
    }

    const std::filesystem::path directory(timeseries_dir);
    std::filesystem::create_directories(directory);

    const std::filesystem::path output_path = directory / (excitation.CaseId + ".csv");

    std::ofstream file(output_path);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open timeseries CSV: " + output_path.string());
    }

    KellerMiksisModel model(parameters, excitation);

    file << "t_s,R_m,U_m_s,T_g_K,p_g_pa,p_inf_pa,delta_T_m,Q_out_W\n";

    for (const BubbleSample& sample : samples) {
        const BubbleState& state = sample.State;

        const double pg_pa = model.GasPressurePa(state);
        const double pinf_pa = model.ExternalPressurePa(sample.TimeS);
        const double delta_t_m = model.ThermalLayerThicknessM(state);
        const double q_out_w = model.HeatLossRateW(state);

        file << sample.TimeS << ','
             << state.R << ','
             << state.U << ','
             << state.Tg << ','
             << pg_pa << ','
             << pinf_pa << ','
             << delta_t_m << ','
             << q_out_w
             << '\n';
    }
}

void WriteTimeSeriesCsvFiles(const std::string& timeseries_dir, const std::vector<BubbleExcitationInput>& inputs, const BubblePhysicalParameters& parameters,
                            const std::vector<CandidateSimulationOutput>& outputs) {
    const std::size_t count = std::min(inputs.size(), outputs.size());

    for (std::size_t i = 0; i < count; ++i) {
        if (outputs[i].Result.Status != "NumericallyInvalid" && !outputs[i].Samples.empty()) {
            WriteBubbleTimeSeriesCsv(timeseries_dir, inputs[i], parameters, outputs[i].Samples);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    try {
        if (argc < 5) {
            std::cerr
                << "Usage: ./solve_bubble_candidates "
                << "<config.yaml> <bubble_excitations.csv> "
                << "<bubble_results_all.csv> <timeseries_dir> [max_cases]\n";
            return 1;
        }

        const std::string config_path = argv[1];
        const std::string input_csv_path = argv[2];
        const std::string output_csv_path = argv[3];
        const std::string timeseries_dir = argv[4];

        const ProjectConfig project_config = YamlReader::ReadProjectConfig(config_path);
        const BubblePhysicalParameters physical_parameters = MakePhysicalParametersFromConfig(project_config);

        std::size_t max_cases = 0;
        if (argc >= 6) {
            max_cases = static_cast<std::size_t>(std::stoull(argv[5]));
        }

        const std::vector<BubbleExcitationInput> inputs = ReadBubbleExcitationsCsv(input_csv_path, max_cases);

        std::cout << "Loaded bubble excitations: " << inputs.size() << '\n';

        for (const BubbleExcitationInput& input : inputs) {
            std::cout
                << input.CaseId
                << " omega=" << input.OmegaRadS
                << " amplitude=" << input.DriveAmplitudePa
                << " phase=" << input.DrivePhaseRad
                << " R0=" << input.EquilibriumRadiusM
                << '\n';
        }

        std::vector<BubbleSimulationResult> simulation_results;
        std::vector<CandidateSimulationOutput> simulation_outputs;

        simulation_results.reserve(inputs.size());
        simulation_outputs.reserve(inputs.size());

        for (const BubbleExcitationInput& input : inputs) {
            CandidateSimulationOutput output = SimulateBubbleCandidate(input, physical_parameters);

            simulation_results.push_back(output.Result);

            std::cout
                << "result "
                << output.Result.CaseId
                << " status=" << output.Result.Status
                << " R_min=" << output.Result.RMinM
                << " R_max=" << output.Result.RMaxM
                << " T_max=" << output.Result.TMaxK
                << " K_R=" << output.Result.CompressionRatio
                << " sbs_candidate=" << output.Result.ThermoMechanicalLuminescenceCandidate
                << '\n';

            simulation_outputs.push_back(std::move(output));
        }

        WriteBubbleResultsCsv(output_csv_path, simulation_results);
        std::cout << "Wrote bubble results to " << output_csv_path << '\n';

        const std::string best_candidates_path = MakeBestCandidatesPath(output_csv_path);
        WriteBestSbsCandidatesCsv(best_candidates_path, simulation_results, 10);
        std::cout << "Wrote best SBS candidates to " << best_candidates_path << '\n';

        WriteTimeSeriesCsvFiles(timeseries_dir, inputs, physical_parameters, simulation_outputs);
        std::cout << "Wrote time series to " << timeseries_dir << '\n';

        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Fatal error: " << exception.what() << '\n';
        return 1;
    }
}