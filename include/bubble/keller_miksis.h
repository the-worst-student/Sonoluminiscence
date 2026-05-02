#pragma once

#include "bubble/bubble_simulation_result.hpp"
#include "io/bubble_excitation_reader.hpp"

struct BubblePhysicalParameters {
    double LiquidDensityKgPerM3 = 1000.0;
    double LiquidSoundSpeedMPerS = 1500.0;
    double LiquidDynamicViscosityPaS = 0.001;
    double LiquidSurfaceTensionNPerM = 0.072;

    double GasMolarConstantJPerMolK = 8.314462618;
    double GasThermalConductivityWPerMK = 0.0177;
    double GasThermalDiffusivityM2PerS = 2.0e-5;
    double GasIsochoricHeatCapacityJPerMolK = 12.471693927;
    double ThermalLayerCoefficient = 2.0;

    double HardCoreRatio = 0.25;
};

class KellerMiksisModel {
public:
    KellerMiksisModel(
        const BubblePhysicalParameters& parameters,
        const BubbleExcitationInput& excitation
    );

    BubbleStateDerivative Evaluate(double t, const BubbleState& state) const;

    double GasPressurePa(const BubbleState& state) const;
    double ExternalPressurePa(double t) const;
    double ExternalPressureDerivativePaS(double t) const;

    double ThermalLayerThicknessM(const BubbleState& state) const;
    double HeatLossRateW(const BubbleState& state) const;

    double HardCoreRadiusM() const;
    double ArgonAmountMol() const;

private:
    double EffectiveVolumeM3(double radius_m) const;
    double EffectiveVolumeDerivativeM3S(const BubbleState& state) const;

    BubblePhysicalParameters parameters_;
    BubbleExcitationInput excitation_;

    double equilibrium_radius_m_ = 0.0;
    double hard_core_radius_m_ = 0.0;
    double argon_amount_mol_ = 0.0;
};