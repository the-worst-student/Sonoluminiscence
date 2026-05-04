#pragma once

#include <complex>
#include <string>

struct BubbleExcitation {
    std::string CaseId = "base_case";

    double FrequencyHz = 0.0;
    double AngularFrequencyRadS = 0.0;

    double BubbleRM = 0.0;
    double BubbleZM = 0.0;

    double StaticPressurePa = 0.0;
    double LiquidTemperatureK = 0.0;
    double EquilibriumRadiusM = 0.0;

    std::complex<double> BubblePressurePa = {0.0, 0.0};

    double PressureAmplitudePa() const;
    double PressureArgumentRad() const;
    double PressurePhaseRad() const;

    static std::string CsvHeader();
    std::string CsvRow() const;
};