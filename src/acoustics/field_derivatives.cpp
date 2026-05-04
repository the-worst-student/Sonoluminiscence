#include <acoustics/field_derivatives.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace {

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
        throw std::runtime_error("Degenerate triangle while computing pressure gradient");
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

std::complex<double> PressureAtNode(
    const AcousticFieldData& field_data,
    const int node_index) {
    const AcousticNodalField& Field =
        field_data.NodalFields.at(static_cast<std::size_t>(node_index));
    return {Field.PressureRealPa, Field.PressureImagPa};
}

void AddWeightedGradientToNode(
    AcousticNodalField* field,
    const std::complex<double>& grad_r,
    const std::complex<double>& grad_z,
    const double weight) {
    field->GradPressureRRealPaPerM += weight * std::real(grad_r);
    field->GradPressureRImagPaPerM += weight * std::imag(grad_r);
    field->GradPressureZRealPaPerM += weight * std::real(grad_z);
    field->GradPressureZImagPaPerM += weight * std::imag(grad_z);
}

void UpdateCellPressureGradient(
    AcousticCellField* field,
    const std::complex<double>& grad_r,
    const std::complex<double>& grad_z) {
    field->GradPressureAbsPaPerM =
        std::sqrt(std::norm(grad_r) + std::norm(grad_z));
}

void FinalizeNodalPressureGradient(
    AcousticNodalField* field,
    const double weight) {
    if (weight <= 0.0) {
        return;
    }

    field->GradPressureRRealPaPerM /= weight;
    field->GradPressureRImagPaPerM /= weight;
    field->GradPressureZRealPaPerM /= weight;
    field->GradPressureZImagPaPerM /= weight;

    const std::complex<double> GradR(
        field->GradPressureRRealPaPerM,
        field->GradPressureRImagPaPerM);
    const std::complex<double> GradZ(
        field->GradPressureZRealPaPerM,
        field->GradPressureZImagPaPerM);

    field->GradPressureAbsPaPerM =
        std::sqrt(std::norm(GradR) + std::norm(GradZ));
}

}  // namespace

void FieldDerivatives::AddPressureGradient(AcousticFieldData* field_data) {
    if (field_data == nullptr) {
        throw std::runtime_error("Cannot add pressure gradient to null field data");
    }

    if (field_data->Nodes.size() != field_data->NodalFields.size()) {
        throw std::runtime_error("Node count does not match nodal acoustic field count");
    }

    if (field_data->Elements.size() != field_data->CellFields.size()) {
        throw std::runtime_error("Element count does not match cell acoustic field count");
    }

    std::vector<double> NodalWeights(field_data->Nodes.size(), 0.0);

    for (std::size_t i = 0; i < field_data->Elements.size(); ++i) {
        const TriangleElement& Element = field_data->Elements[i];
        const MeshNode& Node0 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node0));
        const MeshNode& Node1 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node1));
        const MeshNode& Node2 = field_data->Nodes.at(static_cast<std::size_t>(Element.Node2));

        const ShapeGradients Gradients = ComputeShapeGradients(Node0, Node1, Node2);

        const std::complex<double> Pressure0 = PressureAtNode(*field_data, Element.Node0);
        const std::complex<double> Pressure1 = PressureAtNode(*field_data, Element.Node1);
        const std::complex<double> Pressure2 = PressureAtNode(*field_data, Element.Node2);

        const std::complex<double> GradR =
            Pressure0 * Gradients.Grad0.Dr +
            Pressure1 * Gradients.Grad1.Dr +
            Pressure2 * Gradients.Grad2.Dr;

        const std::complex<double> GradZ =
            Pressure0 * Gradients.Grad0.Dz +
            Pressure1 * Gradients.Grad1.Dz +
            Pressure2 * Gradients.Grad2.Dz;

        UpdateCellPressureGradient(&field_data->CellFields[i], GradR, GradZ);

        const double Weight = Gradients.Area;
        const int NodeIds[3] = {Element.Node0, Element.Node1, Element.Node2};
        for (const int NodeId : NodeIds) {
            const std::size_t NodeIndex = static_cast<std::size_t>(NodeId);
            AddWeightedGradientToNode(
                &field_data->NodalFields.at(NodeIndex),
                GradR,
                GradZ,
                Weight);
            NodalWeights.at(NodeIndex) += Weight;
        }
    }

    for (std::size_t i = 0; i < field_data->NodalFields.size(); ++i) {
        FinalizeNodalPressureGradient(&field_data->NodalFields[i], NodalWeights[i]);
    }
}
