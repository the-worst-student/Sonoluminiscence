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
}