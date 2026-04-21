#include <io/csv_writer.hpp>

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

void CsvWriter::WriteBubbleSamples(
    const std::string& output_path,
    const std::vector<BubbleSample>& samples) {
    std::ofstream Output(output_path);
    if (!Output.is_open()) {
        throw std::runtime_error("Failed to open CSV output file: " + output_path);
    }

    Output << std::scientific
           << std::setprecision(std::numeric_limits<double>::max_digits10);

    Output << "time_s,"
           << "radius_m,"
           << "radius_velocity_m_per_s,"
           << "gas_temperature_k,"
           << "external_pressure_pa,"
           << "gas_pressure_pa\n";

    if (!Output) {
        throw std::runtime_error("Failed to write CSV header to file: " + output_path);
    }

    for (const BubbleSample& Sample : samples) {
        Output << Sample.TimeS << ','
               << Sample.State.RadiusM << ','
               << Sample.State.RadiusVelocityMPerS << ','
               << Sample.State.GasTemperatureK << ','
               << Sample.ExternalPressurePa << ','
               << Sample.GasPressurePa << '\n';

        if (!Output) {
            throw std::runtime_error(
                "Failed while writing bubble samples to CSV file: " + output_path);
        }
    }

    Output.close();
    if (!Output) {
        throw std::runtime_error("Failed to finalize CSV output file: " + output_path);
    }
}