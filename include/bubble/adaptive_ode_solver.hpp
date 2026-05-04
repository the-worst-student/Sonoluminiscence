#pragma once

#include "bubble/bubble_simulation_result.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

struct AdaptiveOdeOptions {
    double InitialStepS = 1.0e-9;
    double MinStepS = 1.0e-14;
    double MaxStepS = 1.0e-7;

    double RelativeTolerance = 1.0e-6;
    double AbsoluteToleranceR = 1.0e-12;
    double AbsoluteToleranceU = 1.0e-3;
    double AbsoluteToleranceT = 1.0e-2;

    std::size_t MaxSteps = 1000000;
};

struct AdaptiveOdeResult {
    bool Success = false;
    std::string FailureReason;
    std::vector<BubbleSample> Samples;
};

class AdaptiveOdeSolver {
public:
    AdaptiveOdeResult Solve(const BubbleState& initial_state, double t_start, double t_end, const AdaptiveOdeOptions& options,
                            const std::function<BubbleStateDerivative(double, const BubbleState&)>& rhs) const;
};