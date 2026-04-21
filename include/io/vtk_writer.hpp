#pragma once

#include <string>

#include <acoustics/field_sampler.hpp>

class VtkWriter {
public:
    static void WriteAxisymmetricSlice(
        const std::string& output_path,
        const AcousticFieldData& field_data,
        double bubble_r,
        double bubble_z);

    static void WritePseudoVolume(
        const std::string& output_path,
        const AcousticFieldData& field_data,
        double bubble_r,
        double bubble_z,
        int sectors);
};
