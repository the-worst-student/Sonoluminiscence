#pragma once

<<<<<<< HEAD
#include "bubble/bubble_rhs.hpp"

#include <cstddef>
#include <optional>
#include <vector>

enum class OdeStopReason {
    Completed,
    RadiusNonPos,
    RadiusBelowMin,
    TemperatureNonPos,
    NonFiniteState,
    NumericalFailure,
    CollapseDetected
};

struct OdeSolverOptions {
    double TimeStepS = 1.0e-10;
    double FinalTimeS = 1.0e-5;
    std::size_t OutputStride = 1;
    bool SaveInitialState = true;
    bool StopAtFirstRebound = false;
    double MinRadiusM = 0.0;
    double MinTemperatureK = 1.0;
};


struct BubbleMetrics {
    bool CollapseDetected = false;
    bool PotentialLuminescence = false;

    double InitialRadiusM = 0.0;

    double MinRadiusM = 0.0;
    double MinRadiusTimeS = 0.0;

    double CollapseTimeS = 0.0;
    double TemperatureAtCollapseK = 0.0;
    double GasPressureAtCollapsePa = 0.0;

    double MaxTemperatureK = 0.0;
    double MaxGasPressurePa = 0.0;
    double MaxInwardSpeedMPerS = 0.0;
};

struct OdeSolution {
    std::vector<BubbleSample> Samples;
    BubbleMetrics Metrics;
    OdeStopReason StopReason = OdeStopReason::Completed;
    bool Success = false;
=======
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
>>>>>>> 7de201c (Update project version)
};

class OdeSolver {
public:
<<<<<<< HEAD
    explicit OdeSolver(const BubbleRhs& rhs, const OdeSolverOptions& options)
    OdeSolution Solve(const BubbleState& initial_state) const;

private:
    BubbleState StepRk4(double t, double dt, const BubbleState& state) const;
    BubbleSample BuildSample(double t, const BubbleState& state) const;

    static bool IsFiniteState(const BubbleState& state);


    std::optional<OdeStopReason> CheckPhysicalState(const BubbleState& state) const;

    void UpdateMetrics(const BubbleSample& sample, BubbleMetrics& metrics) const;

    bool DetectCollapse(const BubbleSample& prev,
                        const BubbleSample& curr,
                        BubbleMetrics& metrics) const;

    static BubbleState AddScaled(const BubbleState& state,
                                    const BubbleDerivative& derivative,
                                    double scale);
    static BubbleState CombineRk4(const BubbleState& state,
                                    double dt,
                                    const BubbleDerivative& k1,
                                    const BubbleDerivative& k2,
                                    const BubbleDerivative& k3,
                                    const BubbleDerivative& k4);

    const BubbleRhs& Rhs;
    OdeSolverOptions Options;
=======
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
>>>>>>> 7de201c (Update project version)
};