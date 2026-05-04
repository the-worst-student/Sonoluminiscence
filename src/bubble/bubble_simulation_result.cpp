#include "bubble/bubble_simulation_result.hpp"

#include "bubble/keller_miksis.hpp"
#include "io/bubble_excitation_reader.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

void FillCommonFields(BubbleSimulationResult& result, const BubbleExcitationInput& excitation, const KellerMiksisModel& model) {
    result.CaseId = excitation.CaseId;

    result.FrequencyHz = excitation.FrequencyHz;
    result.OmegaRadS = excitation.OmegaRadS;
    result.DriveAmplitudePa = excitation.DriveAmplitudePa;
    result.DrivePhaseRad = excitation.DrivePhaseRad;

    result.StaticPressurePa = excitation.StaticPressurePa;
    result.LiquidTemperatureK = excitation.LiquidTemperatureK;

    result.R0M = excitation.EquilibriumRadiusM;
    result.HardCoreRadiusM = model.HardCoreRadiusM();
}

bool IsFiniteState(const KmBubbleState& state) {
    return std::isfinite(state.R) && std::isfinite(state.U) && std::isfinite(state.Tg);
}

std::string ClassifyStatus(const BubbleSimulationResult& result) {
    if (result.ThermoMechanicalLuminescenceCandidate) {
        return "ThermoMechanicalLuminescenceCandidate";
    }

    if (result.MechanicalCollapseOk && !result.ThermalHeatingOk) {
        return "MechanicalCollapseOnly";
    }

    if (result.CompressionRatio < 3.0 && result.TMaxK < 1000.0) {
        return "SoftCompression";
    }

    if ((result.CompressionRatio >= 3.0 && result.CompressionRatio < 8.0) ||
        (result.TMaxK >= 1000.0 && result.TMaxK < 8000.0)) {
        return "ModerateCollapse";
    }

    return "ModerateCollapse";
}

}  // namespace

BubbleSimulationResult AnalyzeBubbleSimulationResult(const BubbleExcitationInput& excitation,
                                                     const KellerMiksisModel& model,
                                                     const BubblePhysicalParameters& parameters,
                                                     const std::vector<KmBubbleSample>& samples,
                                                     bool integration_success,
                                                     const std::string& integration_failure_reason) {
    BubbleSimulationResult result;
    FillCommonFields(result, excitation, model);

    if (!integration_success) {
        result.Status = "NumericallyInvalid";
        result.FailureReason = integration_failure_reason;
        return result;
    }

    if (samples.empty()) {
        result.Status = "NumericallyInvalid";
        result.FailureReason = "EmptyTimeSeries";
        return result;
    }

    double r_min = std::numeric_limits<double>::max();
    double r_max = 0.0;
    double t_max = 0.0;
    double pg_max = 0.0;
    double u_max = 0.0;

    bool min_radius_detected = false;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const KmBubbleState& state = samples[i].State;

        if (!IsFiniteState(state)) {
            result.Status = "NumericallyInvalid";
            result.FailureReason = "NonFiniteStateInTimeSeries";
            return result;
        }

        if (state.R <= 1.001 * model.HardCoreRadiusM() || state.Tg <= 0.0) {
            result.Status = "NumericallyInvalid";
            result.FailureReason = "NonPhysicalStateInTimeSeries";
            return result;
        }

        if (std::abs(state.U) > 0.8 * parameters.LiquidSoundSpeedMPerS ||
            std::abs(1.0 - state.U / parameters.LiquidSoundSpeedMPerS) < 0.05) {
            result.Status = "NumericallyInvalid";
            result.FailureReason = "MachSafetyLimitExceededInTimeSeries";
            return result;
        }

        const double pg = model.GasPressurePa(state);

        if (!std::isfinite(pg) || pg <= 0.0) {
            result.Status = "NumericallyInvalid";
            result.FailureReason = "InvalidGasPressureInTimeSeries";
            return result;
        }

        r_min = std::min(r_min, state.R);
        r_max = std::max(r_max, state.R);
        t_max = std::max(t_max, state.Tg);
        pg_max = std::max(pg_max, pg);
        u_max = std::max(u_max, std::abs(state.U));

        if (i > 0) {
            const KmBubbleState& prev = samples[i - 1].State;
            const KmBubbleState& curr = samples[i].State;

            if (prev.U < 0.0 && curr.U >= 0.0 && curr.U != prev.U) {
                const double alpha = -prev.U / (curr.U - prev.U);
                const double interpolated_radius = prev.R + alpha * (curr.R - prev.R);
                r_min = std::min(r_min, interpolated_radius);
                min_radius_detected = true;
            }
        }
    }

    result.RMinM = r_min;
    result.RMaxM = r_max;
    result.TMaxK = t_max;
    result.PgMaxPa = pg_max;
    result.UMaxMS = u_max;

    if (result.RMinM <= 0.0 || result.RMaxM <= 0.0 || result.R0M <= 0.0 || result.HardCoreRadiusM <= 0.0) {
        result.Status = "NumericallyInvalid";
        result.FailureReason = "InvalidMetricRadius";
        return result;
    }

    result.CompressionRatio = result.R0M / result.RMinM;
    result.ExpansionRatio = result.RMaxM / result.R0M;
    result.DynamicRange = result.RMaxM / result.RMinM;
    result.LiquidMach = result.UMaxMS / parameters.LiquidSoundSpeedMPerS;
    result.HardCoreRatio = result.RMinM / result.HardCoreRadiusM;

    result.MinRadiusDetected = min_radius_detected;

    result.MechanicalCollapseOk =
        result.CompressionRatio >= 8.0 &&
        result.ExpansionRatio >= 5.0 &&
        result.DynamicRange >= 40.0 &&
        result.LiquidMach >= 0.03 &&
        result.HardCoreRatio >= 1.2;

    result.ThermalHeatingOk = result.TMaxK >= 8000.0 && result.PgMaxPa >= 1.0e8;

    result.ThermoMechanicalLuminescenceCandidate = result.MinRadiusDetected && result.MechanicalCollapseOk && result.ThermalHeatingOk;

    result.Status = ClassifyStatus(result);
    return result;
}