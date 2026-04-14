#pragma once
#include <core/config.hpp>
#include <coupling/pressure_coupling.hpp>

struct BubbleState {
    double RadiusM = 0.0;
    double RadiusVelocityMPerS = 0.0;
    double GasTemperatureK = 0.0;
};

struct BubbleDerivative {
    double RadiusRateMPerS = 0.0;
    double RadiusAccelerationMPerS2 = 0.0;
    double GasTemperatureRateKPerS = 0.0;
};

struct BubbleSample {
    double TimeS = 0.0;
    BubbleState State;
    double ExternalPressurePa = 0.0;
    double GasPressurePa = 0.0;
};

struct BubbleModelParameters {
    double LiquidDensityKgPerM3 = 0.0;
    double LiquidViscosityPaS = 0.0;
    double SurfaceTensionNPerM = 0.0;
    double StaticPressurePa = 0.0;
    double LiquidTemperatureK = 0.0;

    double AngularFrequencyRadS = 0.0;

    double EquilibriumRadiusM = 0.0;
    double InitialRadiusM = 0.0;
    double InitialVelocityMPerS = 0.0;
    double InitialTemperatureK = 0.0;

    double HardCoreRadiusM = 0.0;
    double GasCvJPerMolK = 0.0;
    double GasThermalConductivityWPerMK = 0.0;
    double ThermalLayerThicknessM = 0.0;
    double GasMolesMol = 0.0;
};

class BubbleRhs {
public:
    explicit BubbleRhs(const ProjectConfig& config);
    BubbleDerivative Evaluate(double time_s, const BubbleState& state) const;
    double ComputeExternalPressure(double time_s) const;
    double ComputeGasPressure(const BubbleState& state) const;
    BubbleState BuildInitialState() const;
    const BubbleModelParameters& GetParameters() const;
private:
    static BubbleModelParameters BuildParameters(const ProjectConfig& config);
    BubbleModelParameters Parameters;
    PressureCoupling ExternalPressure;
};