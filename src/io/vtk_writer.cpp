#include <io/vtk_writer.hpp>

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
constexpr double cPi = 3.14159265358979323846;
constexpr int cVtkTriangle = 5;
constexpr int cVtkWedge = 13;
constexpr int cVtkVertex = 1;

void EnsureOpen(const std::ofstream& output, const std::string& path) {
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open VTK output file: " + path);
    }
}

std::size_t PointIndex(std::size_t sector, std::size_t node, std::size_t node_count) {
    return sector * node_count + node;
}

std::size_t FindClosestNodeIndex(
    const std::vector<MeshNode>& nodes,
    const double bubble_r,
    const double bubble_z) {
    if (nodes.empty()) {
        throw std::runtime_error("Cannot locate bubble marker in empty node array");
    }

    std::size_t best_index = 0;
    double best_distance2 = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double dr = nodes[i].R - bubble_r;
        const double dz = nodes[i].Z - bubble_z;
        const double distance2 = dr * dr + dz * dz;
        if (distance2 < best_distance2) {
            best_distance2 = distance2;
            best_index = i;
        }
    }
    return best_index;
}

}  // namespace

void VtkWriter::WriteAxisymmetricSlice(
    const std::string& output_path,
    const AcousticFieldData& field_data,
    const double bubble_r,
    const double bubble_z) {
    std::ofstream Output(output_path);
    EnsureOpen(Output, output_path);

    const std::size_t PointCount = field_data.Nodes.size();
    const std::size_t CellCount = field_data.Elements.size();
    const std::size_t BubbleNodeIndex = FindClosestNodeIndex(field_data.Nodes, bubble_r, bubble_z);

    Output << "# vtk DataFile Version 3.0\n";
    Output << "axisymmetric_acoustic_slice\n";
    Output << "ASCII\n";
    Output << "DATASET UNSTRUCTURED_GRID\n";

    Output << "POINTS " << PointCount << " double\n";
    for (const MeshNode& Node : field_data.Nodes) {
        Output << Node.R << ' ' << 0.0 << ' ' << Node.Z << '\n';
    }

    Output << "CELLS " << CellCount << ' ' << CellCount * 4 << "\n";
    for (const TriangleElement& Element : field_data.Elements) {
        Output << 3 << ' '
               << Element.Node0 << ' '
               << Element.Node1 << ' '
               << Element.Node2 << '\n';
    }

    Output << "CELL_TYPES " << CellCount << "\n";
    for (std::size_t i = 0; i < CellCount; ++i) {
        Output << cVtkTriangle << '\n';
    }

    Output << "POINT_DATA " << PointCount << "\n";

    auto WritePointScalar = [&](const std::string& Name, auto Getter) {
        Output << "SCALARS " << Name << " double 1\n";
        Output << "LOOKUP_TABLE default\n";
        for (const AcousticNodalField& Field : field_data.NodalFields) {
            Output << Getter(Field) << '\n';
        }
    };

    WritePointScalar("pressure_re_pa", [](const AcousticNodalField& Field) { return Field.PressureRealPa; });
    WritePointScalar("pressure_im_pa", [](const AcousticNodalField& Field) { return Field.PressureImagPa; });
    WritePointScalar("pressure_abs_pa", [](const AcousticNodalField& Field) { return Field.PressureAbsPa; });
    WritePointScalar("pressure_phase_rad", [](const AcousticNodalField& Field) { return Field.PressurePhaseRad; });

    Output << "SCALARS bubble_marker double 1\n";
    Output << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < PointCount; ++i) {
        Output << (i == BubbleNodeIndex ? 1.0 : 0.0) << '\n';
    }

    Output << "CELL_DATA " << CellCount << "\n";
    Output << "SCALARS pressure_abs_cell_pa double 1\n";
    Output << "LOOKUP_TABLE default\n";
    for (const AcousticCellField& Field : field_data.CellFields) {
        Output << Field.PressureAbsPa << '\n';
    }
}

void VtkWriter::WritePseudoVolume(
    const std::string& output_path,
    const AcousticFieldData& field_data,
    const double bubble_r,
    const double bubble_z,
    const int sectors) {
    if (sectors < 3) {
        throw std::runtime_error("Pseudo-volume export requires at least 3 sectors");
    }

    std::ofstream Output(output_path);
    EnsureOpen(Output, output_path);

    const std::size_t NodeCount2D = field_data.Nodes.size();
    const std::size_t CellCount2D = field_data.Elements.size();
    const std::size_t BubbleNodeIndex2D = FindClosestNodeIndex(field_data.Nodes, bubble_r, bubble_z);
    const std::size_t PointCount = NodeCount2D * static_cast<std::size_t>(sectors) + 1;
    const std::size_t CellCount = CellCount2D * static_cast<std::size_t>(sectors) + 1;

    Output << "# vtk DataFile Version 3.0\n";
    Output << "axisymmetric_acoustic_pseudovolume\n";
    Output << "ASCII\n";
    Output << "DATASET UNSTRUCTURED_GRID\n";

    Output << "POINTS " << PointCount << " double\n";
    for (int Sector = 0; Sector < sectors; ++Sector) {
        const double Phi = 2.0 * cPi * static_cast<double>(Sector) / static_cast<double>(sectors);
        const double CosPhi = std::cos(Phi);
        const double SinPhi = std::sin(Phi);

        for (const MeshNode& Node : field_data.Nodes) {
            const double X = Node.R * CosPhi;
            const double Y = Node.R * SinPhi;
            const double Z = Node.Z;
            Output << X << ' ' << Y << ' ' << Z << '\n';
        }
    }
    const std::size_t BubblePointIndex = PointCount - 1;
    Output << bubble_r << ' ' << 0.0 << ' ' << bubble_z << '\n';

    Output << "CELLS " << CellCount << ' ' << (CellCount - 1) * 7 + 2 << "\n";
    for (int Sector = 0; Sector < sectors; ++Sector) {
        const int NextSector = (Sector + 1) % sectors;
        for (const TriangleElement& Element : field_data.Elements) {
            const std::size_t P0 = PointIndex(static_cast<std::size_t>(Sector), static_cast<std::size_t>(Element.Node0), NodeCount2D);
            const std::size_t P1 = PointIndex(static_cast<std::size_t>(Sector), static_cast<std::size_t>(Element.Node1), NodeCount2D);
            const std::size_t P2 = PointIndex(static_cast<std::size_t>(Sector), static_cast<std::size_t>(Element.Node2), NodeCount2D);
            const std::size_t P3 = PointIndex(static_cast<std::size_t>(NextSector), static_cast<std::size_t>(Element.Node0), NodeCount2D);
            const std::size_t P4 = PointIndex(static_cast<std::size_t>(NextSector), static_cast<std::size_t>(Element.Node1), NodeCount2D);
            const std::size_t P5 = PointIndex(static_cast<std::size_t>(NextSector), static_cast<std::size_t>(Element.Node2), NodeCount2D);

            Output << 6 << ' '
                   << P0 << ' ' << P1 << ' ' << P2 << ' '
                   << P3 << ' ' << P4 << ' ' << P5 << '\n';
        }
    }
    Output << 1 << ' ' << BubblePointIndex << '\n';

    Output << "CELL_TYPES " << CellCount << "\n";
    for (std::size_t i = 0; i < CellCount - 1; ++i) {
        Output << cVtkWedge << '\n';
    }
    Output << cVtkVertex << '\n';

    Output << "POINT_DATA " << PointCount << "\n";

    auto WritePointScalar = [&](const std::string& Name, auto Getter) {
        Output << "SCALARS " << Name << " double 1\n";
        Output << "LOOKUP_TABLE default\n";
        for (int Sector = 0; Sector < sectors; ++Sector) {
            (void)Sector;
            for (const AcousticNodalField& Field : field_data.NodalFields) {
                Output << Getter(Field) << '\n';
            }
        }
        Output << 0.0 << '\n';
    };

    WritePointScalar("pressure_re_pa", [](const AcousticNodalField& Field) { return Field.PressureRealPa; });
    WritePointScalar("pressure_im_pa", [](const AcousticNodalField& Field) { return Field.PressureImagPa; });
    WritePointScalar("pressure_abs_pa", [](const AcousticNodalField& Field) { return Field.PressureAbsPa; });
    WritePointScalar("pressure_phase_rad", [](const AcousticNodalField& Field) { return Field.PressurePhaseRad; });

    Output << "SCALARS bubble_marker double 1\n";
    Output << "LOOKUP_TABLE default\n";
    for (int Sector = 0; Sector < sectors; ++Sector) {
        for (std::size_t Node = 0; Node < NodeCount2D; ++Node) {
            Output << (Node == BubbleNodeIndex2D ? 1.0 : 0.0) << '\n';
        }
    }
    Output << 1.0 << '\n';

    Output << "CELL_DATA " << CellCount << "\n";
    Output << "SCALARS pressure_abs_cell_pa double 1\n";
    Output << "LOOKUP_TABLE default\n";
    for (int Sector = 0; Sector < sectors; ++Sector) {
        (void)Sector;
        for (const AcousticCellField& Field : field_data.CellFields) {
            Output << Field.PressureAbsPa << '\n';
        }
    }
    Output << 0.0 << '\n';
}
