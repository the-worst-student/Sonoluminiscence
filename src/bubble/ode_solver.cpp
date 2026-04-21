#include <bubble/ode_solver.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace {

constexpr std::size_t cStateSize = 3;

double Clamp01(const double value) {
    return std::max(0.0, std::min(1.0, value));
}

double ComputeZeroCrossingAlpha(
    const double previous_velocity_m_s,
    const double current_velocity_m_s) {
    const double denominator = current_velocity_m_s - previous_velocity_m_s;

    if (std::abs(denominator) < 1e-30) {
        return 0.5;
    }

    const double alpha = -previous_velocity_m_s / denominator;
    return Clamp01(alpha);
}

double Lerp(
    const double a,
    const double b,
    const double alpha) {
    return a + (b - a) * alpha;
}

BubbleSample BuildSample(
    const BubbleRhs& rhs,
    const double time_s,
    const BubbleState& state) {
    BubbleSample sample;
    sample.TimeS = time_s;
    sample.State = state;
    sample.ExternalPressurePa = rhs.ComputeExternalPressure(time_s);
    sample.GasPressurePa = rhs.ComputeGasPressure(state);
    return sample;
}

} // namespace

OdeSolver::OdeSolver(const ProjectConfig& config)
    : Config(config), Rhs(config) {
}

OdeSolver::OdeSolver(
    const ProjectConfig& config,
    const PressureCoupling& external_pressure)
    : Config(config), Rhs(config, external_pressure) {
}

void OdeSolver::Solve() {
    ValidateIntegrationConfig();

    Result = OdeSolveResult{};

    gsl_set_error_handler_off();

    const BubbleState InitialState = Rhs.BuildInitialState();

    BubbleStopReason ValidationReason = BubbleStopReason::NotStarted;
    if (!TryValidateState(InitialState, &ValidationReason)) {
        Result.Success = false;
        Result.StopReason = ValidationReason;
        return;
    }

    SaveSample(0.0, InitialState);

    const double FinalTimeS = Config.bubble.integration.final_time_s;
    const double BaseTimeStepS = Config.bubble.integration.time_step_s;
    const int OutputEvery = Config.bubble.integration.output_every;

    if (FinalTimeS == 0.0) {
        Result.Success = true;
        Result.StopReason = BubbleStopReason::CompletedFinalTime;
        return;
    }

    double StateValues[cStateSize];
    StateToArray(InitialState, StateValues);

    gsl_odeiv2_system System;
    System.function = &OdeSolver::EvaluateRhs;
    System.jacobian = nullptr;
    System.dimension = cStateSize;
    System.params = this;

    using StepPtr = std::unique_ptr<gsl_odeiv2_step, void (*)(gsl_odeiv2_step*)>;

    StepPtr Step(
        gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, cStateSize),
        &gsl_odeiv2_step_free);

    if (!Step) {
        throw std::runtime_error("Failed to allocate GSL RK4 stepper");
    }

    double TimeS = 0.0;
    int StepIndex = 0;

    double ErrorValues[cStateSize] = {0.0, 0.0, 0.0};
    double DerivativeIn[cStateSize] = {0.0, 0.0, 0.0};
    double DerivativeOut[cStateSize] = {0.0, 0.0, 0.0};

    BubbleSample PreviousStepSample = Result.Samples.back();

    auto AppendSampleIfNew = [this](const BubbleSample& sample) {
        if (Result.Samples.empty() || Result.Samples.back().TimeS < sample.TimeS) {
            Result.Samples.push_back(sample);
        }
    };

    while (TimeS < FinalTimeS) {
        const double CurrentStepS = std::min(BaseTimeStepS, FinalTimeS - TimeS);

        const int RhsStatus = EvaluateRhs(TimeS, StateValues, DerivativeIn, this);
        if (RhsStatus != GSL_SUCCESS) {
            Result.Success = false;
            Result.StopReason = BubbleStopReason::RhsEvaluationFailed;
            return;
        }

        const int Status = gsl_odeiv2_step_apply(
            Step.get(),
            TimeS,
            CurrentStepS,
            StateValues,
            ErrorValues,
            DerivativeIn,
            DerivativeOut,
            &System);

        if (Status != GSL_SUCCESS) {
            Result.Success = false;
            Result.StopReason = BubbleStopReason::GslStepFailed;
            return;
        }

        TimeS += CurrentStepS;
        ++StepIndex;

        const BubbleState CurrentState = ArrayToState(StateValues);

        if (!TryValidateState(CurrentState, &ValidationReason)) {
            Result.Success = false;
            Result.StopReason = ValidationReason;
            return;
        }

        const BubbleSample CurrentStepSample = BuildSample(Rhs, TimeS, CurrentState);

        if (TryDetectCollapse(PreviousStepSample, CurrentStepSample)) {
            AppendSampleIfNew(PreviousStepSample);
            AppendSampleIfNew(CurrentStepSample);

            Result.Success = true;
            Result.StopReason = BubbleStopReason::CollapseDetected;
            return;
        }

        if (StepIndex % OutputEvery == 0 || TimeS >= FinalTimeS) {
            AppendSampleIfNew(CurrentStepSample);
        }

        PreviousStepSample = CurrentStepSample;
    }

    Result.Success = true;
    Result.StopReason = BubbleStopReason::CompletedFinalTime;
}

bool OdeSolver::IsSolved() const {
    return Result.Success;
}

const std::vector<BubbleSample>& OdeSolver::GetSamples() const {
    return Result.Samples;
}

const BubbleSample& OdeSolver::GetFinalSample() const {
    if (Result.Samples.empty()) {
        throw std::runtime_error("Bubble trajectory is empty");
    }

    return Result.Samples.back();
}

const OdeSolveResult& OdeSolver::GetResult() const {
    return Result;
}

const char* OdeSolver::StopReasonToString(const BubbleStopReason reason) {
    switch (reason) {
        case BubbleStopReason::NotStarted:
            return "NotStarted";
        case BubbleStopReason::CompletedFinalTime:
            return "CompletedFinalTime";
        case BubbleStopReason::CollapseDetected:
            return "CollapseDetected";
        case BubbleStopReason::NonFiniteState:
            return "NonFiniteState";
        case BubbleStopReason::RadiusNonPositive:
            return "RadiusNonPositive";
        case BubbleStopReason::RadiusBelowHardCore:
            return "RadiusBelowHardCore";
        case BubbleStopReason::TemperatureNonPositive:
            return "TemperatureNonPositive";
        case BubbleStopReason::RhsEvaluationFailed:
            return "RhsEvaluationFailed";
        case BubbleStopReason::GslStepFailed:
            return "GslStepFailed";
        default:
            return "Unknown";
    }
}

int OdeSolver::EvaluateRhs(
    const double time_s,
    const double state_values[],
    double derivative_values[],
    void* params) {
    if (params == nullptr) {
        return GSL_EFAULT;
    }

    try {
        const auto* Solver = static_cast<const OdeSolver*>(params);
        const BubbleState State = ArrayToState(state_values);
        const BubbleDerivative Derivative = Solver->Rhs.Evaluate(time_s, State);

        derivative_values[0] = Derivative.RadiusRateMPerS;
        derivative_values[1] = Derivative.RadiusAccelerationMPerS2;
        derivative_values[2] = Derivative.GasTemperatureRateKPerS;

        return GSL_SUCCESS;
    } catch (...) {
        return GSL_EFAILED;
    }
}

BubbleState OdeSolver::ArrayToState(const double state_values[]) {
    BubbleState State;
    State.RadiusM = state_values[0];
    State.RadiusVelocityMPerS = state_values[1];
    State.GasTemperatureK = state_values[2];
    return State;
}

void OdeSolver::StateToArray(
    const BubbleState& state,
    double state_values[]) {
    state_values[0] = state.RadiusM;
    state_values[1] = state.RadiusVelocityMPerS;
    state_values[2] = state.GasTemperatureK;
}

void OdeSolver::ValidateIntegrationConfig() const {
    if (Config.bubble.integration.time_step_s <= 0.0) {
        throw std::runtime_error(
            "Bubble integration time_step_s must be positive");
    }

    if (Config.bubble.integration.final_time_s < 0.0) {
        throw std::runtime_error(
            "Bubble integration final_time_s must be non-negative");
    }

    if (Config.bubble.integration.output_every <= 0) {
        throw std::runtime_error(
            "Bubble integration output_every must be positive");
    }
}

bool OdeSolver::TryValidateState(
    const BubbleState& state,
    BubbleStopReason* reason) const {
    if (reason == nullptr) {
        throw std::runtime_error("Bubble stop reason pointer must not be null");
    }

    if (!std::isfinite(state.RadiusM) ||
        !std::isfinite(state.RadiusVelocityMPerS) ||
        !std::isfinite(state.GasTemperatureK)) {
        *reason = BubbleStopReason::NonFiniteState;
        return false;
    }

    if (state.RadiusM <= 0.0) {
        *reason = BubbleStopReason::RadiusNonPositive;
        return false;
    }

    const BubbleModelParameters& Parameters = Rhs.GetParameters();

    if (state.RadiusM <= Parameters.HardCoreRadiusM) {
        *reason = BubbleStopReason::RadiusBelowHardCore;
        return false;
    }

    if (state.GasTemperatureK <= 0.0) {
        *reason = BubbleStopReason::TemperatureNonPositive;
        return false;
    }

    *reason = BubbleStopReason::NotStarted;
    return true;
}

void OdeSolver::SaveSample(
    const double time_s,
    const BubbleState& state) {
    Result.Samples.push_back(BuildSample(Rhs, time_s, state));
}

bool OdeSolver::TryDetectCollapse(
    const BubbleSample& previous_sample,
    const BubbleSample& current_sample) {
    const double PreviousVelocity = previous_sample.State.RadiusVelocityMPerS;
    const double CurrentVelocity = current_sample.State.RadiusVelocityMPerS;

    if (!(PreviousVelocity < 0.0 && CurrentVelocity >= 0.0)) {
        return false;
    }

    const double Alpha = ComputeZeroCrossingAlpha(
        PreviousVelocity,
        CurrentVelocity);

    Result.Collapse.Detected = true;
    Result.Collapse.TimeS = Lerp(
        previous_sample.TimeS,
        current_sample.TimeS,
        Alpha);
    Result.Collapse.RadiusM = Lerp(
        previous_sample.State.RadiusM,
        current_sample.State.RadiusM,
        Alpha);
    Result.Collapse.TemperatureK = Lerp(
        previous_sample.State.GasTemperatureK,
        current_sample.State.GasTemperatureK,
        Alpha);
    Result.Collapse.GasPressurePa = Lerp(
        previous_sample.GasPressurePa,
        current_sample.GasPressurePa,
        Alpha);

    return true;
}