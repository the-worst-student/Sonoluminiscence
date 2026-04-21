#pragma once

#include <vector>

#include <acoustics/helmholtz_solver.hpp>

struct AcousticNodalField {
    double PressureRealPa = 0.0;
    double PressureImagPa = 0.0;
    double PressureAbsPa = 0.0;
    double PressurePhaseRad = 0.0;
};

struct AcousticCellField {
    double PressureAbsPa = 0.0;
};

struct AcousticFieldData {
    std::vector<MeshNode> Nodes;
    std::vector<TriangleElement> Elements;
    std::vector<AcousticNodalField> NodalFields;
    std::vector<AcousticCellField> CellFields;
};

class FieldSampler {
public:
    static AcousticFieldData BuildFromSolver(const HelmholtzSolver& solver);
};
