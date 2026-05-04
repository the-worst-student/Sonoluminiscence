#include "io/bubble_excitation_reader.hpp"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> fields;

    std::stringstream stream(line);
    std::string field;

    while (std::getline(stream, field, ',')) {
        fields.push_back(field);
    }

    return fields;
}

double ParseDouble(const std::string& text, const std::string& column_name) {
    try {
        return std::stod(text);
    } catch (const std::exception& exception) {
        throw std::runtime_error(
            "Failed to parse column '" + column_name + "' as double. Value: '" +
            text + "'. Parser error: " + exception.what()
        );
    }
}

BubbleExcitationInput ParseBubbleExcitationRow(const std::vector<std::string>& fields, std::size_t line_number) {
    constexpr std::size_t kExpectedColumnCount = 14;

    if (fields.size() != kExpectedColumnCount) {
        throw std::runtime_error(
            "Invalid CSV column count at line " + std::to_string(line_number) +
            ". Expected " + std::to_string(kExpectedColumnCount) +
            ", got " + std::to_string(fields.size())
        );
    }

    BubbleExcitationInput input;

    input.CaseId = fields[0];

    input.FrequencyHz = ParseDouble(fields[1], "frequency_hz");
    input.OmegaRadS = ParseDouble(fields[2], "omega_rad_s");

    input.BubbleRM = ParseDouble(fields[3], "bubble_r_m");
    input.BubbleZM = ParseDouble(fields[4], "bubble_z_m");

    input.StaticPressurePa = ParseDouble(fields[5], "static_pressure_pa");
    input.LiquidTemperatureK = ParseDouble(fields[6], "liquid_temperature_k");
    input.EquilibriumRadiusM = ParseDouble(fields[7], "equilibrium_radius_m");

    input.PbRealPa = ParseDouble(fields[8], "pb_real_pa");
    input.PbImagPa = ParseDouble(fields[9], "pb_imag_pa");
    input.PbAbsPa = ParseDouble(fields[10], "pb_abs_pa");
    input.PbArgumentRad = ParseDouble(fields[11], "pb_argument_rad");

    input.DriveAmplitudePa = ParseDouble(fields[12], "drive_amplitude_pa");
    input.DrivePhaseRad = ParseDouble(fields[13], "drive_phase_rad");

    return input;
}

}  // namespace

std::vector<BubbleExcitationInput> ReadBubbleExcitationsCsv(const std::string& file_path, std::size_t max_cases) {
    std::ifstream file(file_path);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open bubble excitations CSV: " + file_path);
    }

    std::string header;
    if (!std::getline(file, header)) {
        throw std::runtime_error("Bubble excitations CSV is empty: " + file_path);
    }

    std::vector<BubbleExcitationInput> inputs;

    std::string line;
    std::size_t line_number = 1;

    while (std::getline(file, line)) {
        ++line_number;

        if (line.empty()) {
            continue;
        }

        const std::vector<std::string> fields = SplitCsvLine(line);
        BubbleExcitationInput input = ParseBubbleExcitationRow(fields, line_number);

        inputs.push_back(input);

        if (max_cases != 0 && inputs.size() >= max_cases) {
            break;
        }
    }

    return inputs;
}