#include "io/bubble_excitation_reader.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <string>
#include <vector>


int main(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    try {
        if (argc < 5) {
            std::cerr
                << "Usage: ./solve_bubble_candidates "
                << "<config.yaml> <bubble_excitations.csv> "
                << "<bubble_results_all.csv> <timeseries_dir> [max_cases]\n";

            return 1;
        }

        const std::string config_path = argv[1];
        const std::string input_csv_path = argv[2];
        const std::string output_csv_path = argv[3];
        const std::string timeseries_dir = argv[4];

        std::size_t max_cases = 0;

        if (argc >= 6) {
            max_cases = static_cast<std::size_t>(std::stoull(argv[5]));
        }

        const std::vector<BubbleExcitationInput> inputs = ReadBubbleExcitationsCsv(input_csv_path, max_cases);

        std::cout << "Loaded bubble excitations: " << inputs.size() << '\n';

        for (const BubbleExcitationInput& input : inputs) {
            std::cout
                << input.CaseId
                << " omega=" << input.OmegaRadS
                << " amplitude=" << input.DriveAmplitudePa
                << " phase=" << input.DrivePhaseRad
                << " R0=" << input.EquilibriumRadiusM
                << '\n';
        }

        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Fatal error: " << exception.what() << '\n';
        return 1;
    }
}