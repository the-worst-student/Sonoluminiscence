#include <io/result_writer.hpp>

#include <fstream>
#include <stdexcept>

void ResultWriter::WriteTextFile(
    const std::string& output_path,
    const std::string& content) {
    std::ofstream Output(output_path);
    if (!Output.is_open()) {
        throw std::runtime_error("Failed to open text output file: " + output_path);
    }

    Output << content;
}
