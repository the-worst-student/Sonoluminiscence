#pragma once

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
};

class OdeSolver {
public:
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
};