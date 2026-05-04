#include <coupling/bubble_excitation.hpp>

#include <cmath>
#include <iomanip>
#include <sstream>

double BubbleExcitation::PressureAmplitudePa() const {
    return std::abs(BubblePressurePa);
}

double BubbleExcitation::PressureArgumentRad() const {
    return std::arg(BubblePressurePa);
}

double BubbleExcitation::PressurePhaseRad() const {
    return -std::arg(BubblePressurePa);
}

std::string BubbleExcitation::CsvHeader() {
    return "case_id,"
           "frequency_hz,"
           "omega_rad_s,"
           "bubble_r_m,"
           "bubble_z_m,"
           "static_pressure_pa,"
           "liquid_temperature_k,"
           "equilibrium_radius_m,"
           "pb_real_pa,"
           "pb_imag_pa,"
           "pb_abs_pa,"
           "pb_argument_rad,"
           "drive_amplitude_pa,"
           "drive_phase_rad\n";
}

std::string BubbleExcitation::CsvRow() const {
    std::ostringstream Output;
    Output << std::setprecision(17);

    Output << CaseId << ','
           << FrequencyHz << ','
           << AngularFrequencyRadS << ','
           << BubbleRM << ','
           << BubbleZM << ','
           << StaticPressurePa << ','
           << LiquidTemperatureK << ','
           << EquilibriumRadiusM << ','
           << std::real(BubblePressurePa) << ','
           << std::imag(BubblePressurePa) << ','
           << PressureAmplitudePa() << ','
           << PressureArgumentRad() << ','
           << PressureAmplitudePa() << ','
           << PressurePhaseRad() << '\n';

    return Output.str();
}