<<<<<<< HEAD
#pragma once

#include <string>
#include <vector>

#include "bubble/bubble_rhs.hpp"

class CsvWriter {
 public:
    static void WriteBubbleSamples(const std::string& path, const std::vector<BubbleSample>& samples);
};
=======
#include <postprocess/metrics.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

constexpr double cNearHardCoreThreshold = 1.20;

BubbleMetrics EvaluateTrajectoryMetricsOnly(
    const std::vector<BubbleSample>& samples,
    const ProjectConfig& config) {
    if (samples.empty()) {
        throw std::runtime_error(
            "Cannot compute bubble metrics for empty sample history");
    }

    BubbleMetrics MetricsValue;

    MetricsValue.InitialRadiusM = config.bubble.initial_radius_m;
    MetricsValue.InitialTemperatureK = config.bubble.initial_temperature_k;
    MetricsValue.StaticPressurePa = config.liquid.static_pressure_pa;
    MetricsValue.HardCoreRadiusM = config.bubble.gas.hard_core_radius_m;
    MetricsValue.SampleCount = samples.size();

    const BubbleSample& FirstSample = samples.front();
    const BubbleSample& LastSample = samples.back();

    MetricsValue.MinRadiusM = FirstSample.State.RadiusM;
    MetricsValue.TimeAtMinRadiusS = FirstSample.TimeS;

    MetricsValue.MaxTemperatureK = FirstSample.State.GasTemperatureK;
    MetricsValue.TimeAtMaxTemperatureS = FirstSample.TimeS;

    MetricsValue.MaxGasPressurePa = FirstSample.GasPressurePa;
    MetricsValue.TimeAtMaxGasPressureS = FirstSample.TimeS;

    if (FirstSample.State.RadiusVelocityMPerS < 0.0) {
        MetricsValue.MaxInwardSpeedMPerS =
            -FirstSample.State.RadiusVelocityMPerS;
        MetricsValue.TimeAtMaxInwardSpeedS = FirstSample.TimeS;
    }

    for (const BubbleSample& Sample : samples) {
        if (Sample.State.RadiusM < MetricsValue.MinRadiusM) {
            MetricsValue.MinRadiusM = Sample.State.RadiusM;
            MetricsValue.TimeAtMinRadiusS = Sample.TimeS;
        }

        if (Sample.State.GasTemperatureK > MetricsValue.MaxTemperatureK) {
            MetricsValue.MaxTemperatureK = Sample.State.GasTemperatureK;
            MetricsValue.TimeAtMaxTemperatureS = Sample.TimeS;
        }

        if (Sample.GasPressurePa > MetricsValue.MaxGasPressurePa) {
            MetricsValue.MaxGasPressurePa = Sample.GasPressurePa;
            MetricsValue.TimeAtMaxGasPressureS = Sample.TimeS;
        }

        if (Sample.State.RadiusVelocityMPerS < 0.0) {
            const double InwardSpeedMPerS =
                -Sample.State.RadiusVelocityMPerS;

            if (InwardSpeedMPerS > MetricsValue.MaxInwardSpeedMPerS) {
                MetricsValue.MaxInwardSpeedMPerS = InwardSpeedMPerS;
                MetricsValue.TimeAtMaxInwardSpeedS = Sample.TimeS;
            }
        }
    }

    MetricsValue.FinalRadiusM = LastSample.State.RadiusM;
    MetricsValue.FinalTemperatureK = LastSample.State.GasTemperatureK;
    MetricsValue.FinalGasPressurePa = LastSample.GasPressurePa;
    MetricsValue.FinalExternalPressurePa = LastSample.ExternalPressurePa;
    MetricsValue.FinalTimeS = LastSample.TimeS;

    if (MetricsValue.MinRadiusM > 0.0) {
        MetricsValue.CompressionRatio =
            MetricsValue.InitialRadiusM / MetricsValue.MinRadiusM;
    }

    if (MetricsValue.HardCoreRadiusM > 0.0) {
        MetricsValue.HardCoreProximityRatio =
            MetricsValue.MinRadiusM / MetricsValue.HardCoreRadiusM;
    }

    if (MetricsValue.InitialTemperatureK > 0.0) {
        MetricsValue.TemperatureGain =
            MetricsValue.MaxTemperatureK / MetricsValue.InitialTemperatureK;
    }

    if (MetricsValue.StaticPressurePa > 0.0) {
        MetricsValue.PressureGain =
            MetricsValue.MaxGasPressurePa / MetricsValue.StaticPressurePa;
    }

    MetricsValue.NearHardCore =
        MetricsValue.HardCoreProximityRatio > 0.0 &&
        MetricsValue.HardCoreProximityRatio <= cNearHardCoreThreshold;

    return MetricsValue;
}

} // namespace

BubbleMetrics Metrics::EvaluateBubbleMetrics(
    const OdeSolveResult& solve_result,
    const ProjectConfig& config) {
    if (solve_result.Samples.empty()) {
        throw std::runtime_error(
            "Cannot compute bubble metrics for empty ODE solve result");
    }

    BubbleMetrics MetricsValue =
        EvaluateTrajectoryMetricsOnly(solve_result.Samples, config);

    MetricsValue.Success = solve_result.Success;
    MetricsValue.StopReason = solve_result.StopReason;

    MetricsValue.CollapseDetected = solve_result.Collapse.Detected;
    if (solve_result.Collapse.Detected) {
        MetricsValue.CollapseTimeS = solve_result.Collapse.TimeS;
        MetricsValue.RadiusAtCollapseM = solve_result.Collapse.RadiusM;
        MetricsValue.TemperatureAtCollapseK =
            solve_result.Collapse.TemperatureK;
        MetricsValue.GasPressureAtCollapsePa =
            solve_result.Collapse.GasPressurePa;
    }

    return MetricsValue;
}

BubbleMetrics Metrics::EvaluateBubbleMetrics(
    const std::vector<BubbleSample>& samples,
    const ProjectConfig& config) {
    BubbleMetrics MetricsValue =
        EvaluateTrajectoryMetricsOnly(samples, config);

    // Legacy fallback: when only the trajectory is available, we can still
    // compute extrema and final values, but we do not reliably know the
    // solver stop reason or the interpolated collapse point.
    MetricsValue.Success = true;
    MetricsValue.StopReason = BubbleStopReason::CompletedFinalTime;
    MetricsValue.CollapseDetected = false;

    return MetricsValue;
}

void Metrics::ApplyLuminescenceCriteria(
    BubbleMetrics& metrics,
    const LuminescenceCriteria& criteria) {
    metrics.StrongCompression =
        metrics.CompressionRatio >= criteria.MinCompressionRatio;

    metrics.StrongHeating =
        metrics.MaxTemperatureK >= criteria.MinTemperatureK;

    metrics.StrongGasPressure =
        metrics.MaxGasPressurePa >= criteria.MinGasPressurePa;

    metrics.StrongInwardSpeed =
        metrics.MaxInwardSpeedMPerS >= criteria.MinInwardSpeedMPerS;

    metrics.PotentialLuminescence =
        metrics.StrongCompression &&
        metrics.StrongHeating &&
        metrics.StrongGasPressure &&
        metrics.StrongInwardSpeed;
}
>>>>>>> 7de201c (Update project version)
