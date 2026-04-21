#pragma once

#include <vector>

#include <bubble/bubble_rhs.hpp>
#include <core/config.hpp>
#include <coupling/pressure_coupling.hpp>

enum class BubbleStopReason {
    NotStarted,
    CompletedFinalTime,
    CollapseDetected,
    NonFiniteState,
    RadiusNonPositive,
    RadiusBelowHardCore,
    TemperatureNonPositive,
    RhsEvaluationFailed,
    GslStepFailed
};

struct BubbleCollapseInfo {
    bool Detected = false;
    double TimeS = 0.0;
    double RadiusM = 0.0;
    double TemperatureK = 0.0;
    double GasPressurePa = 0.0;
};

struct OdeSolveResult {
    bool Success = false;
    BubbleStopReason StopReason = BubbleStopReason::NotStarted;
    BubbleCollapseInfo Collapse;
    std::vector<BubbleSample> Samples;
};

class OdeSolver {
public:
    explicit OdeSolver(const ProjectConfig& config);
    OdeSolver(const ProjectConfig& config, const PressureCoupling& external_pressure);

    void Solve();

    bool IsSolved() const;
    const std::vector<BubbleSample>& GetSamples() const;
    const BubbleSample& GetFinalSample() const;
    const OdeSolveResult& GetResult() const;

    static const char* StopReasonToString(BubbleStopReason reason);

private:
    static int EvaluateRhs(
        double time_s,
        const double state_values[],
        double derivative_values[],
        void* params);

    static BubbleState ArrayToState(const double state_values[]);
    static void StateToArray(
        const BubbleState& state,
        double state_values[]);

    void ValidateIntegrationConfig() const;

    bool TryValidateState(
        const BubbleState& state,
        BubbleStopReason* reason) const;

    void SaveSample(
        double time_s,
        const BubbleState& state);

    bool TryDetectCollapse(
        const BubbleSample& previous_sample,
        const BubbleSample& current_sample);

    ProjectConfig Config;
    BubbleRhs Rhs;
    OdeSolveResult Result;
};