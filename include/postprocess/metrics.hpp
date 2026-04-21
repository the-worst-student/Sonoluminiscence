#pragma once

#include <cstddef>
#include <vector>

#include <bubble/ode_solver.hpp>
#include <core/config.hpp>

struct BubbleMetrics {
    double InitialRadiusM = 0.0;
    double InitialTemperatureK = 0.0;
    double StaticPressurePa = 0.0;
    double HardCoreRadiusM = 0.0;

    double MinRadiusM = 0.0;
    double TimeAtMinRadiusS = 0.0;

    double MaxTemperatureK = 0.0;
    double TimeAtMaxTemperatureS = 0.0;

    double MaxGasPressurePa = 0.0;
    double TimeAtMaxGasPressureS = 0.0;

    double MaxInwardSpeedMPerS = 0.0;
    double TimeAtMaxInwardSpeedS = 0.0;

    double FinalRadiusM = 0.0;
    double FinalTemperatureK = 0.0;
    double FinalGasPressurePa = 0.0;
    double FinalExternalPressurePa = 0.0;
    double FinalTimeS = 0.0;

    double CompressionRatio = 0.0;
    double HardCoreProximityRatio = 0.0;
    double TemperatureGain = 0.0;
    double PressureGain = 0.0;

    bool Success = false;
    BubbleStopReason StopReason = BubbleStopReason::NotStarted;

    bool CollapseDetected = false;
    double CollapseTimeS = 0.0;
    double RadiusAtCollapseM = 0.0;
    double TemperatureAtCollapseK = 0.0;
    double GasPressureAtCollapsePa = 0.0;

    bool StrongCompression = false;
    bool StrongHeating = false;
    bool StrongGasPressure = false;
    bool StrongInwardSpeed = false;
    bool NearHardCore = false;
    bool PotentialLuminescence = false;

    std::size_t SampleCount = 0;
};

struct LuminescenceCriteria {
    double MinTemperatureK = 1500.0;
    double MinGasPressurePa = 5.0e6;
    double MinInwardSpeedMPerS = 0.0;
    double MinCompressionRatio = 5.0;
};

class Metrics {
public:
    static BubbleMetrics EvaluateBubbleMetrics(
        const OdeSolveResult& solve_result,
        const ProjectConfig& config);

    static BubbleMetrics EvaluateBubbleMetrics(
        const std::vector<BubbleSample>& samples,
        const ProjectConfig& config);

    static void ApplyLuminescenceCriteria(
        BubbleMetrics& metrics,
        const LuminescenceCriteria& criteria);
};