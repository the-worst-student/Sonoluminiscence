<<<<<<< HEAD
#include "bubble/ode_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

explicit OdeSolver::OdeSolver(const BubbleRhs& rhs, const OdeSolverOptions& options): Rhs(rhs), Options(options) {
    if (!(Options.TimeStepS > 0.0)) {
        throw std::invalid_argument("OdeSolverOptions.TimeStepS must be positive");
    }
    if (!(Options.FinalTimeS > 0.0)) {
        throw std::invalid_argument("OdeSolverOptions.FinalTimeS must be positive");
    }
    if (Options.OutputStride < 1) {
        throw std::invalid_argument("OdeSolverOptions.OutputStride must be >= 1");
    }
    if (!(Options.MinRadiusM >= 0.0)) {
        throw std::invalid_argument("OdeSolverOptions.MinRadiusM must be non-negative");
    }
    if (!(Options.MinTemperatureK > 0.0)) {
        throw std::invalid_argument("OdeSolverOptions.MinTemperatureK must be positive");
    }
}

bool OdeSolver::IsFiniteState(const BubbleState& state) {
    return std::isfinite(state.RadiusM) && std::isfinite(state.RadiusVelocityMPerS) && std::isfinite(state.GasTemperatureK);
}

std::optional<OdeStopReason> OdeSolver::CheckPhysicalState(const BubbleState& state) const {
    if (!(state.RadiusM > 0.0)) {
        return OdeStopReason::RadiusNonPos;
    }
    if (state.RadiusM < Options.MinRadiusM) {
        return OdeStopReason::RadiusBelowMin;
    }
    if (!(state.GasTemperatureK > Options.MinTemperatureK)) {
        return OdeStopReason::TemperatureNonPos;
    }
    return std::nullopt;
}

BubbleState OdeSolver::AddScaled(const BubbleState& state,
                                 const BubbleDerivative& derivative,
                                 double scale) {
    BubbleState result;
    result.RadiusM = state.RadiusM + scale * derivative.RadiusRateMPerS;
    result.RadiusVelocityMPerS =
        state.RadiusVelocityMPerS + scale * derivative.RadiusAccelerationMPerS2;
    result.GasTemperatureK =
        state.GasTemperatureK + scale * derivative.GasTemperatureRateKPerS;
    return result;
}

BubbleState OdeSolver::CombineRk4(const BubbleState& state,
                                  double dt,
                                  const BubbleDerivative& k1,
                                  const BubbleDerivative& k2,
                                  const BubbleDerivative& k3,
                                  const BubbleDerivative& k4) {
    const double sixth = dt / 6.0;

    BubbleState result;
    result.RadiusM =
        state.RadiusM + sixth * (k1.RadiusRateMPerS
                               + 2.0 * k2.RadiusRateMPerS
                               + 2.0 * k3.RadiusRateMPerS
                               + k4.RadiusRateMPerS);

    result.RadiusVelocityMPerS =
        state.RadiusVelocityMPerS + sixth * (k1.RadiusAccelerationMPerS2
                                           + 2.0 * k2.RadiusAccelerationMPerS2
                                           + 2.0 * k3.RadiusAccelerationMPerS2
                                           + k4.RadiusAccelerationMPerS2);

    result.GasTemperatureK =
        state.GasTemperatureK + sixth * (k1.GasTemperatureRateKPerS
                                       + 2.0 * k2.GasTemperatureRateKPerS
                                       + 2.0 * k3.GasTemperatureRateKPerS
                                       + k4.GasTemperatureRateKPerS);

    return result;
}

BubbleState OdeSolver::StepRk4(double t, double dt, const BubbleState& state) const {
    const BubbleDerivative k1 = Rhs.Evaluate(t, state);
    const BubbleDerivative k2 = Rhs.Evaluate(t + 0.5 * dt, AddScaled(state, k1, 0.5 * dt));
    const BubbleDerivative k3 = Rhs.Evaluate(t + 0.5 * dt, AddScaled(state, k2, 0.5 * dt));
    const BubbleDerivative k4 = Rhs.Evaluate(t + dt, AddScaled(state, k3, dt));

    return CombineRk4(state, dt, k1, k2, k3, k4);
}

BubbleSample OdeSolver::BuildSample(double t, const BubbleState& state) const {
    BubbleSample sample;
    sample.TimeS = t;
    sample.State = state;
    sample.ExternalPressurePa = Rhs.ComputeExternalPressure(t);
    sample.GasPressurePa = Rhs.ComputeGasPressure(state);
    return sample;
}

void OdeSolver::UpdateMetrics(const BubbleSample& sample, BubbleMetrics& metrics) const {
    if (sample.State.RadiusM < metrics.MinRadiusM) {
        metrics.MinRadiusM = sample.State.RadiusM;
        metrics.MinRadiusTimeS = sample.TimeS;
    }

    if (sample.State.GasTemperatureK > metrics.MaxTemperatureK) {
        metrics.MaxTemperatureK = sample.State.GasTemperatureK;
    }

    if (sample.GasPressurePa > metrics.MaxGasPressurePa) {
        metrics.MaxGasPressurePa = sample.GasPressurePa;
    }

    const double inward_speed = -sample.State.RadiusVelocityMPerS;
    if (inward_speed > metrics.MaxInwardSpeedMPerS) {
        metrics.MaxInwardSpeedMPerS = inward_speed;
    }
}

bool OdeSolver::DetectCollapse(const BubbleSample& prev,
                               const BubbleSample& curr,
                               BubbleMetrics& metrics) const {
    if (metrics.CollapseDetected) {
        return false;
    }

    const double v_prev = prev.State.RadiusVelocityMPerS;
    const double v_curr = curr.State.RadiusVelocityMPerS;

    if (!(v_prev < 0.0 && v_curr >= 0.0)) {
        return false;
    }

    const double dv = v_curr - v_prev;
    double alpha = (dv != 0.0) ? (-v_prev / dv) : 0.0;
    alpha = std::max(0.0, std::min(1.0, alpha));

    const double dt = curr.TimeS - prev.TimeS;
    const double t_collapse = prev.TimeS + alpha * dt;

    const double temperature_at_collapse =
        prev.State.GasTemperatureK
        + alpha * (curr.State.GasTemperatureK - prev.State.GasTemperatureK);

    const double pressure_at_collapse =
        prev.GasPressurePa + alpha * (curr.GasPressurePa - prev.GasPressurePa);

    const double radius_at_collapse =
        prev.State.RadiusM + 0.5 * v_prev * alpha * dt;

    if (radius_at_collapse < metrics.MinRadiusM) {
        metrics.MinRadiusM = radius_at_collapse;
        metrics.MinRadiusTimeS = t_collapse;
    }

    metrics.CollapseDetected = true;
    metrics.CollapseTimeS = t_collapse;
    metrics.TemperatureAtCollapseK = temperature_at_collapse;
    metrics.GasPressureAtCollapsePa = pressure_at_collapse;

    return true;
}

OdeSolution OdeSolver::Solve(const BubbleState& initial_state) const {
    OdeSolution solution;
    solution.StopReason = OdeStopReason::Completed;
    solution.Success = false;

    if (!IsFiniteState(initial_state)) {
        solution.StopReason = OdeStopReason::NonFiniteState;
        return solution;
    }

    if (auto bad = CheckPhysicalState(initial_state)) {
        solution.StopReason = *bad;
        return solution;
    }

    BubbleMetrics& metrics = solution.Metrics;
    metrics.InitialRadiusM = initial_state.RadiusM;
    metrics.MinRadiusM = initial_state.RadiusM;
    metrics.MinRadiusTimeS = 0.0;
    metrics.MaxTemperatureK = initial_state.GasTemperatureK;
    metrics.MaxGasPressurePa = 0.0;
    metrics.MaxInwardSpeedMPerS =
        std::max(0.0, -initial_state.RadiusVelocityMPerS);

    BubbleSample prev_sample;
    try {
        prev_sample = BuildSample(0.0, initial_state);
    } catch (...) {
        solution.StopReason = OdeStopReason::NumericalFailure;
        return solution;
    }

    metrics.MaxGasPressurePa = prev_sample.GasPressurePa;

    if (Options.SaveInitialState) {
        solution.Samples.push_back(prev_sample);
    }

    double t = 0.0;
    BubbleState state = initial_state;
    std::size_t step_index = 0;
    bool last_sample_stored = Options.SaveInitialState;

    while (t < Options.FinalTimeS) {
        const double dt = std::min(Options.TimeStepS, Options.FinalTimeS - t);
        const double t_next = t + dt;
        ++step_index;

        BubbleSample curr_sample;
        try {
            const BubbleState next_state = StepRk4(t, dt, state);

            if (!IsFiniteState(next_state)) {
                solution.StopReason = OdeStopReason::NonFiniteState;
                return solution;
            }

            if (auto bad = CheckPhysicalState(next_state)) {
                solution.StopReason = *bad;
                return solution;
            }

            curr_sample = BuildSample(t_next, next_state);
            state = next_state;
        } catch (...) {
            solution.StopReason = OdeStopReason::NumericalFailure;
            return solution;
        }

        UpdateMetrics(curr_sample, metrics);
        const bool collapse_now = DetectCollapse(prev_sample, curr_sample, metrics);

        const bool is_final_step = !(t_next < Options.FinalTimeS);
        const bool on_stride = (step_index % Options.OutputStride == 0);

        if (on_stride || is_final_step) {
            solution.Samples.push_back(curr_sample);
            last_sample_stored = true;
        } else {
            last_sample_stored = false;
        }

        prev_sample = curr_sample;
        t = t_next;

        if (collapse_now && Options.StopAtFirstRebound) {
            if (!last_sample_stored) {
                solution.Samples.push_back(curr_sample);
            }
            solution.StopReason = OdeStopReason::CollapseDetected;
            solution.Success = true;
            return solution;
        }
    }

    solution.StopReason = OdeStopReason::Completed;
    solution.Success = true;
    return solution;
=======
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
>>>>>>> 7de201c (Update project version)
}