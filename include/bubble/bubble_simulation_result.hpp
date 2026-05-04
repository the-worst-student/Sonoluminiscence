#pragma once

#include <string>
#include <vector>

struct BubbleState {
    double R = 0.0;
    double U = 0.0;
    double Tg = 0.0;
};

struct BubbleStateDerivative {
    double dRdt = 0.0;
    double dUdt = 0.0;
    double dTgdt = 0.0;
};

struct BubbleSample {
    double TimeS = 0.0;
    BubbleState State;
};

struct BubbleSimulationResult {
    std::string CaseId;
    std::string Status;
    std::string FailureReason;

    double FrequencyHz = 0.0;
    double OmegaRadS = 0.0;
    double DriveAmplitudePa = 0.0;
    double DrivePhaseRad = 0.0;

    double StaticPressurePa = 0.0;
    double LiquidTemperatureK = 0.0;

    double R0M = 0.0;
    double HardCoreRadiusM = 0.0;

    double RMinM = 0.0;
    double RMaxM = 0.0;

    double TMaxK = 0.0;
    double PgMaxPa = 0.0;
    double UMaxMS = 0.0;

    double CompressionRatio = 0.0;
    double ExpansionRatio = 0.0;
    double DynamicRange = 0.0;
    double LiquidMach = 0.0;
    double HardCoreRatio = 0.0;

    bool MinRadiusDetected = false;
    bool MechanicalCollapseOk = false;
    bool ThermalHeatingOk = false;
    bool ThermoMechanicalLuminescenceCandidate = false;
};

struct BubbleExcitationInput;
struct BubblePhysicalParameters;
class KellerMiksisModel;

BubbleSimulationResult AnalyzeBubbleSimulationResult(const BubbleExcitationInput& excitation,
                                                     const KellerMiksisModel& model,
                                                     const BubblePhysicalParameters& parameters,
                                                     const std::vector<BubbleSample>& samples,
                                                     bool integration_success,
                                                     const std::string& integration_failure_reason);