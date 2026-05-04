#include "bubble/adaptive_ode_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <string>

namespace {

BubbleState AddScaledState(const BubbleState& state, const BubbleStateDerivative& derivative, double scale) {
    BubbleState result;
    result.R = state.R + scale * derivative.dRdt;
    result.U = state.U + scale * derivative.dUdt;
    result.Tg = state.Tg + scale * derivative.dTgdt;
    return result;
}

BubbleState AddWeightedDerivatives(const BubbleState& state,
                                   const BubbleStateDerivative& k1, double b1,
                                   const BubbleStateDerivative& k2, double b2,
                                   const BubbleStateDerivative& k3, double b3,
                                   const BubbleStateDerivative& k4, double b4,
                                   const BubbleStateDerivative& k5, double b5,
                                   const BubbleStateDerivative& k6, double b6,
                                   double dt) {
    BubbleState result;

    result.R = state.R + dt * (b1 * k1.dRdt + b2 * k2.dRdt + b3 * k3.dRdt +
                               b4 * k4.dRdt + b5 * k5.dRdt + b6 * k6.dRdt);

    result.U = state.U + dt * (b1 * k1.dUdt + b2 * k2.dUdt + b3 * k3.dUdt +
                               b4 * k4.dUdt + b5 * k5.dUdt + b6 * k6.dUdt);

    result.Tg = state.Tg + dt * (b1 * k1.dTgdt + b2 * k2.dTgdt + b3 * k3.dTgdt +
                                 b4 * k4.dTgdt + b5 * k5.dTgdt + b6 * k6.dTgdt);

    return result;
}

struct Rk45StepResult {
    BubbleState HighOrderState;
    BubbleState LowOrderState;
    double ErrorNorm = 0.0;
};

double ErrorScale(double absolute_tolerance, double relative_tolerance, double old_value, double new_value) {
    const double magnitude = std::max(std::abs(old_value), std::abs(new_value));
    return absolute_tolerance + relative_tolerance * magnitude;
}

double EstimateErrorNorm(const BubbleState& previous_state,
                         const BubbleState& high_state,
                         const BubbleState& low_state,
                         const AdaptiveOdeOptions& options) {
    const double r_scale = ErrorScale(options.AbsoluteToleranceR, options.RelativeTolerance, previous_state.R, high_state.R);
    const double u_scale = ErrorScale(options.AbsoluteToleranceU, options.RelativeTolerance, previous_state.U, high_state.U);
    const double tg_scale = ErrorScale(options.AbsoluteToleranceT, options.RelativeTolerance, previous_state.Tg, high_state.Tg);

    const double r_error = std::abs(high_state.R - low_state.R) / r_scale;
    const double u_error = std::abs(high_state.U - low_state.U) / u_scale;
    const double tg_error = std::abs(high_state.Tg - low_state.Tg) / tg_scale;

    return std::max(r_error, std::max(u_error, tg_error));
}

Rk45StepResult MakeRk45Step(const BubbleState& state, double t, double dt,
                            const AdaptiveOdeOptions& options,
                            const std::function<BubbleStateDerivative(double, const BubbleState&)>& rhs) {
    const BubbleStateDerivative k1 = rhs(t, state);

    const BubbleState y2 = AddScaledState(state, k1, dt * 1.0 / 5.0);
    const BubbleStateDerivative k2 = rhs(t + dt * 1.0 / 5.0, y2);

    const BubbleState y3 = AddWeightedDerivatives(state, k1, 3.0 / 40.0, k2, 9.0 / 40.0,
                                                  k1, 0.0, k1, 0.0, k1, 0.0, k1, 0.0, dt);
    const BubbleStateDerivative k3 = rhs(t + dt * 3.0 / 10.0, y3);

    const BubbleState y4 = AddWeightedDerivatives(state, k1, 3.0 / 10.0, k2, -9.0 / 10.0,
                                                  k3, 6.0 / 5.0, k1, 0.0, k1, 0.0, k1, 0.0, dt);
    const BubbleStateDerivative k4 = rhs(t + dt * 3.0 / 5.0, y4);

    const BubbleState y5 = AddWeightedDerivatives(state, k1, -11.0 / 54.0, k2, 5.0 / 2.0,
                                                  k3, -70.0 / 27.0, k4, 35.0 / 27.0,
                                                  k1, 0.0, k1, 0.0, dt);
    const BubbleStateDerivative k5 = rhs(t + dt, y5);

    const BubbleState y6 = AddWeightedDerivatives(state, k1, 1631.0 / 55296.0, k2, 175.0 / 512.0,
                                                  k3, 575.0 / 13824.0, k4, 44275.0 / 110592.0,
                                                  k5, 253.0 / 4096.0, k1, 0.0, dt);
    const BubbleStateDerivative k6 = rhs(t + dt * 7.0 / 8.0, y6);

    const BubbleState high_state = AddWeightedDerivatives(state, k1, 37.0 / 378.0, k2, 0.0,
                                                          k3, 250.0 / 621.0, k4, 125.0 / 594.0,
                                                          k5, 0.0, k6, 512.0 / 1771.0, dt);

    const BubbleState low_state = AddWeightedDerivatives(state, k1, 2825.0 / 27648.0, k2, 0.0,
                                                         k3, 18575.0 / 48384.0, k4, 13525.0 / 55296.0,
                                                         k5, 277.0 / 14336.0, k6, 1.0 / 4.0, dt);

    Rk45StepResult result;
    result.HighOrderState = high_state;
    result.LowOrderState = low_state;
    result.ErrorNorm = EstimateErrorNorm(state, high_state, low_state, options);
    return result;
}

bool IsFiniteState(const BubbleState& state) {
    return std::isfinite(state.R) && std::isfinite(state.U) && std::isfinite(state.Tg);
}

double ComputeStepFactor(double error_norm) {
    constexpr double kSafety = 0.9;
    constexpr double kMinFactor = 0.2;
    constexpr double kMaxFactor = 5.0;

    if (error_norm <= 1.0e-16) {
        return kMaxFactor;
    }

    const double factor = kSafety * std::pow(error_norm, -0.2);
    return std::clamp(factor, kMinFactor, kMaxFactor);
}

}  // namespace

AdaptiveOdeResult AdaptiveOdeSolver::Solve(const BubbleState& initial_state, double t_start, double t_end, const AdaptiveOdeOptions& options,
                                           const std::function<BubbleStateDerivative(double, const BubbleState&)>& rhs) const {
    AdaptiveOdeResult result;

    if (t_end <= t_start) {
        result.FailureReason = "InvalidTimeInterval";
        return result;
    }

    if (options.InitialStepS <= 0.0 || options.MinStepS <= 0.0 || options.MaxStepS <= 0.0) {
        result.FailureReason = "InvalidStepOptions";
        return result;
    }

    if (options.MinStepS > options.MaxStepS) {
        result.FailureReason = "MinStepGreaterThanMaxStep";
        return result;
    }

    if (options.RelativeTolerance <= 0.0 || options.AbsoluteToleranceR <= 0.0 ||
        options.AbsoluteToleranceU <= 0.0 || options.AbsoluteToleranceT <= 0.0) {
        result.FailureReason = "InvalidToleranceOptions";
        return result;
    }

    if (options.MaxSteps == 0) {
        result.FailureReason = "InvalidMaxSteps";
        return result;
    }

    if (!IsFiniteState(initial_state) || initial_state.R <= 0.0 || initial_state.Tg <= 0.0) {
        result.FailureReason = "InvalidInitialState";
        return result;
    }

    double t = t_start;
    double dt = std::clamp(options.InitialStepS, options.MinStepS, options.MaxStepS);
    BubbleState state = initial_state;

    result.Samples.push_back(BubbleSample{t, state});

    for (std::size_t step_index = 0; step_index < options.MaxSteps && t < t_end; ++step_index) {
        const bool final_step = (t + dt >= t_end);

        if (final_step) {
            dt = t_end - t;
        }

        if (dt <= 0.0) {
            result.FailureReason = "NonPositiveStep";
            return result;
        }

        if (dt < options.MinStepS && !final_step) {
            result.FailureReason = "StepUnderflow";
            return result;
        }

        Rk45StepResult step_result;

        try {
            step_result = MakeRk45Step(state, t, dt, options, rhs);
        } catch (const std::exception& exception) {
            result.FailureReason = std::string("RhsEvaluationFailed: ") + exception.what();
            return result;
        }

        if (!std::isfinite(step_result.ErrorNorm)) {
            result.FailureReason = "NonFiniteErrorNorm";
            return result;
        }

        if (!IsFiniteState(step_result.HighOrderState)) {
            result.FailureReason = "NonFiniteState";
            return result;
        }

        if (step_result.HighOrderState.R <= 0.0 || step_result.HighOrderState.Tg <= 0.0) {
            result.FailureReason = "NonPhysicalState";
            return result;
        }

        const double factor = ComputeStepFactor(step_result.ErrorNorm);

        if (step_result.ErrorNorm <= 1.0) {
            t += dt;
            state = step_result.HighOrderState;
            result.Samples.push_back(BubbleSample{t, state});
        }

        dt = std::clamp(dt * factor, options.MinStepS, options.MaxStepS);
    }

    if (t < t_end) {
        result.FailureReason = "MaxStepsExceeded";
        return result;
    }

    result.Success = true;
    return result;
}