#pragma once

#include <complex>
#include <cstddef>
#include <string>

#include <core/config.hpp>
#include <postprocess/metrics.hpp>

struct PipelineOutputs {
    std::string MeshPath = "mesh.msh";
    std::string BubbleCsvPath = "bubble_solution.csv";
    std::string SliceVtkPath = "acoustic_slice.vtk";
    std::string PseudoVolumeVtkPath = "acoustic_pseudovolume.vtk";
    std::string SummaryPath = "pipeline_summary.txt";
    std::string BubbleFramesDirectory = "bubble_frames";
    int PseudoVolumeSectors = 48;
};

struct PipelineResult {
    std::complex<double> BubblePressurePa = {0.0, 0.0};
    double BubblePressureAbsPa = 0.0;
    double BubblePressurePhaseRad = 0.0;
    std::size_t SavedSamples = 0;
    BubbleMetrics Metrics;
    PipelineOutputs Outputs;
};

class Pipeline {
public:
    explicit Pipeline(const ProjectConfig& config);

    PipelineResult Run(const PipelineOutputs& outputs) const;

private:
    ProjectConfig Config;
};