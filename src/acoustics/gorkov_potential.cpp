#include <acoustics/gorkov_potential.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace {
constexpr double cPi = 3.14159265358979323846;

struct ShapeGradients {
    Gradient2D Grad0;
    Gradient2D Grad1;
    Gradient2D Grad2;
    double Area = 0.0;
};

ShapeGradients ComputeShapeGradients(
    const MeshNode& node0,
    const MeshNode& node1,
    const MeshNode& node2) {
    const double Denominator =
        (node1.Z - node2.Z) * (node0.R - node2.R) +
        (node2.R - node1.R) * (node0.Z - node2.Z);

    if (std::abs(Denominator) < 1e-18) {
        throw std::runtime_error("Degenerate triangle while computing Gorkov force");
    }

    ShapeGradients Result;
    Result.Grad0.Dr = (node1.Z - node2.Z) / Denominator;
    Result.Grad0.Dz = (node2.R - node1.R) / Denominator;
    Result.Grad1.Dr = (node2.Z - node0.Z) / Denominator;
    Result.Grad1.Dz = (node0.R - node2.R) / Denominator;
    Result.Grad2.Dr = -Result.Grad0.Dr - Result.Grad1.Dr;
    Result.Grad2.Dz = -Result.Grad0.Dz - Result.Grad1.Dz;
    Result.Area = 0.5 * std::abs(Denominator);
    return Result;
}

void ComputeAcousticVelocity(
    AcousticNodalField* field,
    const GorkovParameters& parameters) {
    const double Denominator =
        parameters.AngularFrequencyRadS * parameters.LiquidDensityKgM3;
    if (Denominator <= 0.0) {
        throw std::runtime_error("Invalid omega*rho while computing acoustic velocity");
    }

    const std::complex<double> InvIOverOmegaRho(0.0, -1.0 / Denominator);
    const std::complex<double> GradR(
        field->GradPressureRRealPaPerM,
        field->GradPressureRImagPaPerM);
    const std::complex<double> GradZ(
        field->GradPressureZRealPaPerM,
        field->GradPressureZImagPaPerM);

    const std::complex<double> VelocityR = GradR * InvIOverOmegaRho;
    const std::complex<double> VelocityZ = GradZ * InvIOverOmegaRho;

    field->VelocityRRealMS = std::real(VelocityR);
    field->VelocityRImagMS = std::imag(VelocityR);
    field->VelocityZRealMS = std::real(VelocityZ);
    field->VelocityZImagMS = std::imag(VelocityZ);
    field->VelocityAbsMS = std::sqrt(std::norm(VelocityR) + std::norm(VelocityZ));
}

void ComputeGorkovPotentialAtNode(
    AcousticNodalField* field,
    const GorkovParameters& parameters) {
    const double BubbleVolume =
        4.0 * cPi * parameters.BubbleRadiusM * parameters.BubbleRadiusM *
        parameters.BubbleRadiusM / 3.0;

    const double PressureEnergyTerm =
        parameters.CompressibilityContrast * field->PressureAbsPa * field->PressureAbsPa /
        (4.0 * parameters.LiquidDensityKgM3 *
         parameters.SoundSpeedMS * parameters.SoundSpeedMS);

    const double KineticEnergyTerm =
        parameters.DensityContrast * 3.0 * parameters.LiquidDensityKgM3 *
        field->VelocityAbsMS * field->VelocityAbsMS / 8.0;

    field->GorkovPotentialJ = BubbleVolume * (PressureEnergyTerm - KineticEnergyTerm);
}

void AddWeightedForceToNode(
    AcousticNodalField* field,
    const double force_r,
    const double force_z,
    const double weight) {
    field->GorkovForceRN += weight * force_r;
    field->GorkovForceZN += weight * force_z;
}

void UpdateCellAverages(
    AcousticFieldData* field_data) {
    for (std::size_t i = 0; i < field_data->Elements.size(); ++i) {
        const TriangleElement& Element = field_data->Elements[i];
        const AcousticNodalField& Field0 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node0));
        const AcousticNodalField& Field1 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node1));
        const AcousticNodalField& Field2 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node2));

        AcousticCellField& Cell = field_data->CellFields[i];
        Cell.VelocityAbsMS =
            (Field0.VelocityAbsMS + Field1.VelocityAbsMS + Field2.VelocityAbsMS) / 3.0;
        Cell.GorkovPotentialJ =
            (Field0.GorkovPotentialJ + Field1.GorkovPotentialJ + Field2.GorkovPotentialJ) / 3.0;
        Cell.GorkovForceAbsN =
            (Field0.GorkovForceAbsN + Field1.GorkovForceAbsN + Field2.GorkovForceAbsN) / 3.0;
        Cell.CandidateScore =
            (Field0.CandidateScore + Field1.CandidateScore + Field2.CandidateScore) / 3.0;
    }
}

}  // namespace

GorkovParameters GorkovPotential::FromConfig(
    const ProjectConfig& config,
    const double angular_frequency_rad_s) {
    GorkovParameters Parameters;
    Parameters.BubbleRadiusM = config.bubble.equilibrium_radius_m;
    Parameters.LiquidDensityKgM3 = config.liquid.density_kg_m3;
    Parameters.SoundSpeedMS = config.liquid.sound_speed_m_s;
    Parameters.AngularFrequencyRadS = angular_frequency_rad_s;
    Parameters.CompressibilityContrast = 1.0;
    Parameters.DensityContrast = -2.0;
    return Parameters;
}

void GorkovPotential::AddGorkovFields(
    AcousticFieldData* field_data,
    const GorkovParameters& parameters) {
    if (field_data == nullptr) {
        throw std::runtime_error("Cannot add Gorkov fields to null field data");
    }

    if (parameters.BubbleRadiusM <= 0.0 ||
        parameters.LiquidDensityKgM3 <= 0.0 ||
        parameters.SoundSpeedMS <= 0.0 ||
        parameters.AngularFrequencyRadS <= 0.0) {
        throw std::runtime_error("Invalid parameters for Gorkov potential computation");
    }

    for (AcousticNodalField& Field : field_data->NodalFields) {
        ComputeAcousticVelocity(&Field, parameters);
        ComputeGorkovPotentialAtNode(&Field, parameters);
    }

    std::vector<double> NodalWeights(field_data->Nodes.size(), 0.0);

    for (std::size_t i = 0; i < field_data->Elements.size(); ++i) {
        const TriangleElement& Element = field_data->Elements[i];
        const MeshNode& Node0 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node0));
        const MeshNode& Node1 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node1));
        const MeshNode& Node2 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node2));
        const ShapeGradients Gradients = ComputeShapeGradients(Node0, Node1, Node2);

        const double Potential0 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node0)).GorkovPotentialJ;
        const double Potential1 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node1)).GorkovPotentialJ;
        const double Potential2 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node2)).GorkovPotentialJ;

        const double GradUR =
            Potential0 * Gradients.Grad0.Dr +
            Potential1 * Gradients.Grad1.Dr +
            Potential2 * Gradients.Grad2.Dr;
        const double GradUZ =
            Potential0 * Gradients.Grad0.Dz +
            Potential1 * Gradients.Grad1.Dz +
            Potential2 * Gradients.Grad2.Dz;

        const double ForceR = -GradUR;
        const double ForceZ = -GradUZ;

        const int NodeIds[3] = {Element.Node0, Element.Node1, Element.Node2};
        for (const int NodeId : NodeIds) {
            const std::size_t NodeIndex = static_cast<std::size_t>(NodeId);
            AddWeightedForceToNode(
                &field_data->NodalFields.at(NodeIndex),
                ForceR,
                ForceZ,
                Gradients.Area);
            NodalWeights.at(NodeIndex) += Gradients.Area;
        }
    }

    for (std::size_t i = 0; i < field_data->NodalFields.size(); ++i) {
        AcousticNodalField& Field = field_data->NodalFields[i];
        const double Weight = NodalWeights[i];
        if (Weight > 0.0) {
            Field.GorkovForceRN /= Weight;
            Field.GorkovForceZN /= Weight;
        }
        Field.GorkovForceAbsN =
            std::sqrt(Field.GorkovForceRN * Field.GorkovForceRN +
                      Field.GorkovForceZN * Field.GorkovForceZN);
    }

    UpdateCellAverages(field_data);
}
