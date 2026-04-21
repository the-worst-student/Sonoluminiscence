#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include <coupling/pipeline.hpp>
#include <io/yaml_reader.hpp>

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

}  // namespace

int main(const int argc, char** argv) {
    const std::string ConfigPath =
        argc > 1 ? ResolveExistingPath(argv[1])
                 : ResolveExistingPath("configs/base.yaml");

    PipelineOutputs Outputs;
    if (argc > 2) {
        Outputs.MeshPath = argv[2];
    }
    if (argc > 3) {
        Outputs.BubbleCsvPath = argv[3];
    }
    if (argc > 4) {
        Outputs.SliceVtkPath = argv[4];
    }
    if (argc > 5) {
        Outputs.PseudoVolumeVtkPath = argv[5];
    }
    if (argc > 6) {
        Outputs.SummaryPath = argv[6];
    }
    if (argc > 7) {
        Outputs.BubbleFramesDirectory = argv[7];
    }
    if (argc > 8) {
        Outputs.PseudoVolumeSectors = std::stoi(argv[8]);
    }

    try {
        const ProjectConfig Config = YamlReader::ReadProjectConfig(ConfigPath);
        Pipeline FullPipeline(Config);
        const PipelineResult Result = FullPipeline.Run(Outputs);
        const BubbleMetrics& Metrics = Result.Metrics;

        std::cout << "Full E2E pipeline finished\n";
        std::cout << "Config: " << ConfigPath << '\n';
        std::cout << "Mesh: " << Result.Outputs.MeshPath << '\n';
        std::cout << "Bubble CSV: " << Result.Outputs.BubbleCsvPath << '\n';
        std::cout << "Acoustic slice VTK: " << Result.Outputs.SliceVtkPath << '\n';
        std::cout << "Acoustic pseudovolume VTK: "
                  << Result.Outputs.PseudoVolumeVtkPath << '\n';
        std::cout << "Summary: " << Result.Outputs.SummaryPath << '\n';
        std::cout << "Bubble frames: " << Result.Outputs.BubbleFramesDirectory << '\n';

        std::cout << "Pb abs [Pa]: " << Result.BubblePressureAbsPa << '\n';
        std::cout << "Pb phase [rad]: " << Result.BubblePressurePhaseRad << '\n';
        std::cout << "Saved bubble samples: " << Result.SavedSamples << '\n';

        std::cout << "Bubble solve success: "
                  << std::boolalpha << Metrics.Success << '\n';
        std::cout << "Bubble stop reason: "
                  << OdeSolver::StopReasonToString(Metrics.StopReason) << '\n';

        std::cout << "Min radius [m]: " << Metrics.MinRadiusM << '\n';
        std::cout << "Time at min radius [s]: " << Metrics.TimeAtMinRadiusS << '\n';
        std::cout << "Max gas temperature [K]: " << Metrics.MaxTemperatureK << '\n';
        std::cout << "Max gas pressure [Pa]: " << Metrics.MaxGasPressurePa << '\n';
        std::cout << "Max inward speed [m/s]: "
                  << Metrics.MaxInwardSpeedMPerS << '\n';
        std::cout << "Compression ratio R0/Rmin: "
                  << Metrics.CompressionRatio << '\n';

        std::cout << "Collapse detected: "
                  << std::boolalpha << Metrics.CollapseDetected << '\n';
        if (Metrics.CollapseDetected) {
            std::cout << "Collapse time [s]: " << Metrics.CollapseTimeS << '\n';
            std::cout << "Radius at collapse [m]: "
                      << Metrics.RadiusAtCollapseM << '\n';
            std::cout << "Temperature at collapse [K]: "
                      << Metrics.TemperatureAtCollapseK << '\n';
            std::cout << "Gas pressure at collapse [Pa]: "
                      << Metrics.GasPressureAtCollapsePa << '\n';
        }

        std::cout << "Potential luminescence: "
                  << std::boolalpha << Metrics.PotentialLuminescence << '\n';

        return Metrics.Success ? 0 : 2;
    } catch (const std::exception& Error) {
        std::cerr << "Error while running full pipeline: "
                  << Error.what() << '\n';
        return 1;
    }
}