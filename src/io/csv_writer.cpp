#include "io/csv_writer.hpp"

#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <stdexcept>

void CsvWriter::WriteBubbleSamples(const std::string& path, const std::vector<BubbleSample>& samples) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error(
            "CsvWriter: cannot open file for writing: " + path);
    }

    file << std::scientific;
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    file << "time_s,radius_m,radius_velocity_m_per_s," "gas_temperature_k,external_pressure_pa,gas_pressure_pa\n";

    for (const BubbleSample& sample : samples) {
        file << sample.TimeS << ','
             << sample.State.RadiusM << ','
             << sample.State.RadiusVelocityMPerS << ','
             << sample.State.GasTemperatureK << ','
             << sample.ExternalPressurePa << ','
             << sample.GasPressurePa << '\n';
    }

    if (!file) {
        throw std::runtime_error("CsvWriter: write failed for file: " + path);
    }
}