 #pragma once

#include <cstddef>
#include <string>
#include <vector>

struct BubbleExcitationInput {
    std::string CaseId;

    double FrequencyHz = 0.0;
    double OmegaRadS = 0.0;

    double BubbleRM = 0.0;
    double BubbleZM = 0.0;

    double StaticPressurePa = 0.0;
    double LiquidTemperatureK = 0.0;
    double EquilibriumRadiusM = 0.0;

    double PbRealPa = 0.0;
    double PbImagPa = 0.0;
    double PbAbsPa = 0.0;
    double PbArgumentRad = 0.0;

    double DriveAmplitudePa = 0.0;
    double DrivePhaseRad = 0.0;
};

std::vector<BubbleExcitationInput> ReadBubbleExcitationsCsv(
    const std::string& file_path,
    std::size_t max_cases = 0
);