#pragma once

#include <vector>

#include <acoustics/helmholtz_solver.hpp>

struct AcousticNodalField {
    double PressureRealPa = 0.0;
    double PressureImagPa = 0.0;
    double PressureAbsPa = 0.0;
    double PressurePhaseRad = 0.0;

    double GradPressureRRealPaPerM = 0.0;
    double GradPressureRImagPaPerM = 0.0;
    double GradPressureZRealPaPerM = 0.0;
    double GradPressureZImagPaPerM = 0.0;
    double GradPressureAbsPaPerM = 0.0;

    double VelocityRRealMS = 0.0;
    double VelocityRImagMS = 0.0;
    double VelocityZRealMS = 0.0;
    double VelocityZImagMS = 0.0;
    double VelocityAbsMS = 0.0;

    double GorkovPotentialJ = 0.0;
    double GorkovForceRN = 0.0;
    double GorkovForceZN = 0.0;
    double GorkovForceAbsN = 0.0;

    double CandidateMarker = 0.0;
    double CandidateScore = 0.0;
};

struct AcousticCellField {
    double PressureAbsPa = 0.0;
    double GradPressureAbsPaPerM = 0.0;
    double VelocityAbsMS = 0.0;
    double GorkovPotentialJ = 0.0;
    double GorkovForceAbsN = 0.0;
    double CandidateScore = 0.0;
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
