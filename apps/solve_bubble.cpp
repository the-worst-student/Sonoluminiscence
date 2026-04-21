#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <bubble/ode_solver.hpp>
#include <core/config.hpp>
#include <io/csv_writer.hpp>
#include <io/yaml_reader.hpp>
#include <postprocess/metrics.hpp>

namespace {

std::string ResolveExistingPath(const std::string& path) {
    namespace fs = std::filesystem;

    const fs::path Direct(path);
    if (fs::exists(Direct)) {
        return Direct.string();
    }

    const fs::path Parent = fs::path("..") / Direct;
    if (fs::exists(Parent)) {
        return Parent.string();
    }

    return path;
}

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