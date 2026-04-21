#include <acoustics/field_sampler.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>

AcousticFieldData FieldSampler::BuildFromSolver(const HelmholtzSolver& solver) {
    if (!solver.IsSolved()) {
        throw std::runtime_error("Cannot sample acoustic field before solve");
    }

    const std::vector<MeshNode>& Nodes = solver.GetNodes();
    const std::vector<TriangleElement>& Elements = solver.GetElements();
    const HelmholtzSolution& Solution = solver.GetSolution();

    if (Nodes.size() != Solution.NodalValues.size()) {
        throw std::runtime_error("Acoustic solution size does not match node count");
    }

    AcousticFieldData Data;
    Data.Nodes = Nodes;
    Data.Elements = Elements;
    Data.NodalFields.resize(Nodes.size());
    Data.CellFields.resize(Elements.size());

    for (std::size_t i = 0; i < Nodes.size(); ++i) {
        const std::complex<double> Value = Solution.NodalValues[i];
        AcousticNodalField Field;
        Field.PressureRealPa = std::real(Value);
        Field.PressureImagPa = std::imag(Value);
        Field.PressureAbsPa = std::abs(Value);
        Field.PressurePhaseRad = std::arg(Value);
        Data.NodalFields[i] = Field;
    }

    for (std::size_t i = 0; i < Elements.size(); ++i) {
        const TriangleElement& Element = Elements[i];
        const double Value0 = Data.NodalFields[static_cast<std::size_t>(Element.Node0)].PressureAbsPa;
        const double Value1 = Data.NodalFields[static_cast<std::size_t>(Element.Node1)].PressureAbsPa;
        const double Value2 = Data.NodalFields[static_cast<std::size_t>(Element.Node2)].PressureAbsPa;

        AcousticCellField Field;
        Field.PressureAbsPa = (Value0 + Value1 + Value2) / 3.0;
        Data.CellFields[i] = Field;
    }

    return Data;
}
