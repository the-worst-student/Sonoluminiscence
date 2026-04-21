#include <exception>
#include <filesystem>
#include <iostream>
<<<<<<< HEAD
#include <stdexcept>
#include <string>

#include "bubble/bubble_rhs.hpp"
#include "bubble/ode_solver.hpp"
#include "io/csv_writer.hpp"
#include "io/yaml_reader.hpp"
#include "postprocess/metrics.hpp"
=======
#include <string>
#include <vector>

#include <bubble/ode_solver.hpp>
#include <core/config.hpp>
#include <io/csv_writer.hpp>
#include <io/yaml_reader.hpp>
#include <postprocess/metrics.hpp>
>>>>>>> 7de201c (Update project version)

namespace {

std::string ResolveExistingPath(const std::string& path) {
    namespace fs = std::filesystem;

<<<<<<< HEAD
    const fs::path direct(path);
    if (fs::exists(direct)) {
        return direct.string();
    }

    const fs::path parent = fs::path("..") / direct;
    if (fs::exists(parent)) {
        return parent.string();
=======
    const fs::path Direct(path);
    if (fs::exists(Direct)) {
        return Direct.string();
    }

    const fs::path Parent = fs::path("..") / Direct;
    if (fs::exists(Parent)) {
        return Parent.string();
>>>>>>> 7de201c (Update project version)
    }

    return path;
}

<<<<<<< HEAD
const char* StopReasonToString(OdeStopReason reason) {
    switch (reason) {
        case OdeStopReason::Completed:           return "Completed";
        case OdeStopReason::RadiusNonPos:        return "RadiusNonPos";
        case OdeStopReason::RadiusBelowMin:      return "RadiusBelowMin";
        case OdeStopReason::TemperatureNonPos:   return "TemperatureNonPos";
        case OdeStopReason::NonFiniteState:      return "NonFiniteState";
        case OdeStopReason::NumericalFailure:    return "NumericalFailure";
        case OdeStopReason::CollapseDetected:    return "CollapseDetected";
    }
    return "Unknown";
}

OdeSolverOptions BuildSolverOptions(const ProjectConfig& config) {
    const BubbleIntegrationConfig& integration = config.bubble.integration;
    if (integration.output_every < 1) {
        throw std::invalid_argument(
            "bubble.integration.output_every must be >= 1");
    }

    OdeSolverOptions options;
    options.TimeStepS = integration.time_step_s;
    options.FinalTimeS = integration.final_time_s;
    options.OutputStride = static_cast<std::size_t>(integration.output_every);
    options.MinRadiusM = config.bubble.gas.hard_core_radius_m;
    return options;
}

}

int main(const int argc, char** argv) {
    const std::string config_path =
        argc > 1 ? ResolveExistingPath(argv[1]) : ResolveExistingPath("configs/base.yaml");
    const std::string output_csv_path =
        argc > 2 ? std::string(argv[2]) : std::string("bubble_trajectory.csv");

    try {
        const ProjectConfig config = YamlReader::ReadProjectConfig(config_path);

        const BubbleRhs rhs(config);
        const OdeSolverOptions options = BuildSolverOptions(config);
        const OdeSolver solver(rhs, options);

        const BubbleState initial_state = rhs.BuildInitialState();

        std::cout << "Start Solve\n";
        OdeSolution solution = solver.Solve(initial_state);
        std::cout << "Solve done\n";

        const LuminescenceCriteria criteria;
        ApplyLuminescenceCriteria(solution.Metrics, criteria);

        const BubbleMetrics& m = solution.Metrics;

        std::cout << "Bubble solver finished\n";
        std::cout << "Config: " << config_path << '\n';
        std::cout << "Stop reason: " << StopReasonToString(solution.StopReason) << '\n';
        std::cout << "Success: " << std::boolalpha << solution.Success << '\n';
        std::cout << "Samples stored: " << solution.Samples.size() << '\n';

        std::cout << "Initial radius [m]: " << m.InitialRadiusM << '\n';
        std::cout << "Min radius [m]: " << m.MinRadiusM << '\n';
        std::cout << "Min radius time [s]: " << m.MinRadiusTimeS << '\n';
        std::cout << "Max temperature [K]: " << m.MaxTemperatureK << '\n';
        std::cout << "Max gas pressure [Pa]: " << m.MaxGasPressurePa << '\n';
        std::cout << "Max inward speed [m/s]: " << m.MaxInwardSpeedMPerS << '\n';

        std::cout << "Collapse detected: " << std::boolalpha << m.CollapseDetected << '\n';
        if (m.CollapseDetected) {
            std::cout << "Collapse time [s]: " << m.CollapseTimeS << '\n';
            std::cout << "Temperature at collapse [K]: " << m.TemperatureAtCollapseK << '\n';
            std::cout << "Gas pressure at collapse [Pa]: " << m.GasPressureAtCollapsePa << '\n';
        }

        std::cout << "Potential luminescence: "
                  << std::boolalpha << m.PotentialLuminescence << '\n';

        CsvWriter::WriteBubbleSamples(output_csv_path, solution.Samples);
        std::cout << "Trajectory written: " << output_csv_path << '\n';

        if (!solution.Success) {
            return 2;
        }
    } catch (const std::exception& error) {
        std::cerr << "Error while solving bubble: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
=======
} // namespace

int main(const int argc, char** argv) {
    const std::string ConfigPath =
        argc > 1 ? ResolveExistingPath(argv[1])
                 : ResolveExistingPath("configs/base.yaml");

    const std::string OutputCsvPath =
        argc > 2 ? argv[2] : "bubble_solution.csv";

    try {
        const ProjectConfig Config = YamlReader::ReadProjectConfig(ConfigPath);

        OdeSolver Solver(Config);
        Solver.Solve();

        const OdeSolveResult& SolveResult = Solver.GetResult();
        const std::vector<BubbleSample>& Samples = SolveResult.Samples;

        if (!Samples.empty()) {
            CsvWriter::WriteBubbleSamples(OutputCsvPath, Samples);
        }

        std::cout << "Bubble solver finished\n";
        std::cout << "Config: " << ConfigPath << '\n';
        std::cout << "CSV output: " << OutputCsvPath << '\n';
        std::cout << "Success: " << (SolveResult.Success ? "true" : "false") << '\n';
        std::cout << "Stop reason: "
                  << OdeSolver::StopReasonToString(SolveResult.StopReason) << '\n';
        std::cout << "Saved samples: " << Samples.size() << '\n';

        if (Samples.empty()) {
            std::cout << "No trajectory samples were produced\n";
            return SolveResult.Success ? 0 : 2;
        }

        BubbleMetrics Summary =
            Metrics::EvaluateBubbleMetrics(SolveResult, Config);

        const LuminescenceCriteria Criteria;
        Metrics::ApplyLuminescenceCriteria(Summary, Criteria);

        const BubbleSample& FinalSample = Samples.back();

        std::cout << "Final time [s]: " << FinalSample.TimeS << '\n';
        std::cout << "Final radius [m]: " << FinalSample.State.RadiusM << '\n';
        std::cout << "Final radius velocity [m/s]: "
                  << FinalSample.State.RadiusVelocityMPerS << '\n';
        std::cout << "Final gas temperature [K]: "
                  << FinalSample.State.GasTemperatureK << '\n';
        std::cout << "Final external pressure [Pa]: "
                  << FinalSample.ExternalPressurePa << '\n';
        std::cout << "Final gas pressure [Pa]: "
                  << FinalSample.GasPressurePa << '\n';

        std::cout << "Min radius [m]: " << Summary.MinRadiusM << '\n';
        std::cout << "Time at min radius [s]: " << Summary.TimeAtMinRadiusS << '\n';
        std::cout << "Max gas temperature [K]: " << Summary.MaxTemperatureK << '\n';
        std::cout << "Max gas pressure [Pa]: " << Summary.MaxGasPressurePa << '\n';
        std::cout << "Max inward speed [m/s]: " << Summary.MaxInwardSpeedMPerS << '\n';
        std::cout << "Compression ratio [-]: " << Summary.CompressionRatio << '\n';
        std::cout << "Potential luminescence: "
                  << (Summary.PotentialLuminescence ? "true" : "false") << '\n';

        std::cout << "Collapse detected: "
                  << (Summary.CollapseDetected ? "true" : "false") << '\n';

        if (Summary.CollapseDetected) {
            std::cout << "Collapse time [s]: " << Summary.CollapseTimeS << '\n';
            std::cout << "Radius at collapse [m]: "
                      << Summary.RadiusAtCollapseM << '\n';
            std::cout << "Temperature at collapse [K]: "
                      << Summary.TemperatureAtCollapseK << '\n';
            std::cout << "Gas pressure at collapse [Pa]: "
                      << Summary.GasPressureAtCollapsePa << '\n';
        }

        return SolveResult.Success ? 0 : 2;
    } catch (const std::exception& Error) {
        std::cerr << "Error while solving bubble ODE: "
                  << Error.what() << '\n';
        return 1;
    }
}
>>>>>>> 7de201c (Update project version)
