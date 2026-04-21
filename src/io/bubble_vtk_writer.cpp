#include <io/bubble_vtk_writer.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double cPi = 3.14159265358979323846;

struct Vec3 {
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
};

void EnsureSegments(const BubbleVtkWriter::Options& options) {
    if (options.AzimuthalSegments < 3) {
        throw std::runtime_error("Bubble animation requires at least 3 azimuthal segments");
    }
    if (options.PolarSegments < 3) {
        throw std::runtime_error("Bubble animation requires at least 3 polar segments");
    }
}

std::string FramePath(const std::string& dir, std::size_t index) {
    std::ostringstream name;
    name << dir << "/bubble_frame_" << std::setw(4) << std::setfill('0') << index << ".vtk";
    return name.str();
}

void WriteFrame(
    const std::string& output_path,
    const BubbleSample& sample,
    double bubble_r_m,
    double bubble_z_m,
    int azimuthal_segments,
    int polar_segments) {
    std::ofstream output(output_path);
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open bubble VTK output file: " + output_path);
    }

    const std::size_t point_count = static_cast<std::size_t>((polar_segments + 1) * azimuthal_segments);
    const std::size_t polygon_count = static_cast<std::size_t>(polar_segments * azimuthal_segments);

    std::vector<Vec3> points;
    points.reserve(point_count);

    const double radius = sample.State.RadiusM;
    for (int i = 0; i <= polar_segments; ++i) {
        const double theta = cPi * static_cast<double>(i) / static_cast<double>(polar_segments);
        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        for (int j = 0; j < azimuthal_segments; ++j) {
            const double phi = 2.0 * cPi * static_cast<double>(j) / static_cast<double>(azimuthal_segments);
            const double cos_phi = std::cos(phi);
            const double sin_phi = std::sin(phi);
            Vec3 p;
            p.X = bubble_r_m + radius * sin_theta * cos_phi;
            p.Y = radius * sin_theta * sin_phi;
            p.Z = bubble_z_m + radius * cos_theta;
            points.push_back(p);
        }
    }

    output << "# vtk DataFile Version 3.0\n";
    output << "bubble_animation_frame\n";
    output << "ASCII\n";
    output << "DATASET POLYDATA\n";

    output << "POINTS " << points.size() << " double\n";
    for (const Vec3& p : points) {
        output << p.X << ' ' << p.Y << ' ' << p.Z << '\n';
    }

    output << "POLYGONS " << polygon_count << ' ' << polygon_count * 5 << "\n";
    for (int i = 0; i < polar_segments; ++i) {
        for (int j = 0; j < azimuthal_segments; ++j) {
            const int next_j = (j + 1) % azimuthal_segments;
            const std::size_t p00 = static_cast<std::size_t>(i * azimuthal_segments + j);
            const std::size_t p01 = static_cast<std::size_t>(i * azimuthal_segments + next_j);
            const std::size_t p10 = static_cast<std::size_t>((i + 1) * azimuthal_segments + j);
            const std::size_t p11 = static_cast<std::size_t>((i + 1) * azimuthal_segments + next_j);
            output << 4 << ' ' << p00 << ' ' << p01 << ' ' << p11 << ' ' << p10 << '\n';
        }
    }

    output << "POINT_DATA " << points.size() << "\n";
    output << "SCALARS gas_temperature_k double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.State.GasTemperatureK << '\n';
    }

    output << "SCALARS gas_pressure_pa double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.GasPressurePa << '\n';
    }

    output << "SCALARS external_pressure_pa double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.ExternalPressurePa << '\n';
    }

    output << "SCALARS radius_m double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.State.RadiusM << '\n';
    }

    output << "SCALARS time_s double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.TimeS << '\n';
    }

    output << "SCALARS temperature_colormap_hint double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.State.GasTemperatureK << '\n';
    }
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.GasPressurePa << '\n';
    }

    output << "SCALARS external_pressure_pa double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        output << sample.ExternalPressurePa << '\n';
    }
}

}  // namespace

void BubbleVtkWriter::WriteAnimationFrames(
    const std::vector<BubbleSample>& samples,
    const double bubble_r_m,
    const double bubble_z_m,
    const Options& options) {
    EnsureSegments(options);
    if (samples.empty()) {
        throw std::runtime_error("Cannot write bubble animation for empty trajectory");
    }

    std::filesystem::create_directories(options.OutputDirectory);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        WriteFrame(
            FramePath(options.OutputDirectory, i),
            samples[i],
            bubble_r_m,
            bubble_z_m,
            options.AzimuthalSegments,
            options.PolarSegments);
    }
}
