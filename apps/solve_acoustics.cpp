#include <complex>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "acoustics/acoustics_problem.hpp"
#include "acoustics/boundary_conditions.hpp"
#include "acoustics/helmholtz_solver.hpp"
#include "coupling/pressure_coupling.hpp"
#include "io/yaml_reader.hpp"
#include <acoustics/field_sampler.hpp>
#include <coupling/bubble_excitation.hpp>
#include <io/result_writer.hpp>
#include <io/vtk_writer.hpp>

namespace {

std::string ResolveExistingPath(const std::string& path) {
    namespace fs = std::filesystem;

    const fs::path direct(path);
    if (fs::exists(direct)) {
        return direct.string();
    }

    const fs::path parent = fs::path("..") / direct;
    if (fs::exists(parent)) {
        return parent.string();
    }

    return path;
}

}  // namespace

int main(const int argc, char** argv) {
    const std::string config_path =
        argc > 1 ? ResolveExistingPath(argv[1]) : ResolveExistingPath("configs/base.yaml");
    const std::string mesh_path =
        argc > 2 ? ResolveExistingPath(argv[2]) : ResolveExistingPath("mesh.msh");
    const std::string excitation_path =
        argc > 3 ? argv[3] : "bubble_excitation.csv";
    const std::string slice_vtk_path =
        argc > 4 ? argv[4] : "acoustic_slice.vtk";
    const std::string pseudovolume_vtk_path =
        argc > 5 ? argv[5] : "acoustic_pseudovolume.vtk";
    const int pseudovolume_sectors =
        argc > 6 ? std::stoi(argv[6]) : 48;
    try {
        const ProjectConfig config = YamlReader::ReadProjectConfig(config_path);
        const AcousticsProblem problem = AcousticsProblem::FromConfig(config);
        const BoundaryConditions conditions(problem);

        HelmholtzSolver solver(problem, conditions);

        std::cout << "Start LoadMesh\n";
        solver.LoadMesh(mesh_path);
        std::cout << "LoadMesh done\n";

        std::cout << "Start Assemble\n";
        solver.Assemble();
        std::cout << "Assemble done\n";

        std::cout << "Start Solve\n";
        solver.Solve();
        std::cout << "Solve done\n";

        const double bubble_r = config.geometry.bubble_position.r_m;
        const double bubble_z = config.geometry.bubble_position.z_m;
        const std::complex<double> pb = solver.SamplePressure(bubble_r, bubble_z);

        BubbleExcitation excitation;
        excitation.CaseId =
            config.project.run_id.empty() ? "base_case" : config.project.run_id;
        excitation.FrequencyHz = config.acoustics.frequency_hz;
        excitation.AngularFrequencyRadS = problem.Omega;
        excitation.BubbleRM = bubble_r;
        excitation.BubbleZM = bubble_z;
        excitation.StaticPressurePa = config.liquid.static_pressure_pa;
        excitation.LiquidTemperatureK = config.liquid.temperature_k;
        excitation.EquilibriumRadiusM = config.bubble.equilibrium_radius_m;
        excitation.BubblePressurePa = pb;

        const PressureCoupling coupling(
            excitation.StaticPressurePa,
            excitation.PressureAmplitudePa(),
            excitation.AngularFrequencyRadS,
            excitation.PressurePhaseRad());
        const double period_s = 2.0 * 3.14159265358979323846 /
                                coupling.GetAngularFrequencyRadS();

        const std::string csv_content = BubbleExcitation::CsvHeader() + excitation.CsvRow();
        ResultWriter::WriteTextFile(excitation_path, csv_content);
        const AcousticFieldData field_data = FieldSampler::BuildFromSolver(solver);
        VtkWriter::WriteAxisymmetricSlice(
        slice_vtk_path,
        field_data,
        bubble_r,
        bubble_z);

        VtkWriter::WritePseudoVolume(
        pseudovolume_vtk_path,
        field_data,
        bubble_r,
        bubble_z,
        pseudovolume_sectors);

        std::cout << "Acoustics solver finished\n";
        std::cout << "Config: " << config_path << '\n';
        std::cout << "Mesh: " << mesh_path << '\n';
        std::cout << "Mesh loaded: " << std::boolalpha << solver.IsMeshLoaded() << '\n';
        std::cout << "Assembled: " << std::boolalpha << solver.IsAssembled() << '\n';
        std::cout << "Solved: " << std::boolalpha << solver.IsSolved() << '\n';
        std::cout << "Bubble point r [m]: " << bubble_r << '\n';
        std::cout << "Bubble point z [m]: " << bubble_z << '\n';
        std::cout << "Pb real [Pa]: " << std::real(pb) << '\n';
        std::cout << "Pb imag [Pa]: " << std::imag(pb) << '\n';
        std::cout << "Pb abs [Pa]: " << std::abs(pb) << '\n';
        std::cout << "Pb argument [rad]: " << excitation.PressureArgumentRad() << '\n';
        std::cout << "Drive phase [rad]: "
                  << excitation.PressurePhaseRad() << '\n';
        std::cout << "Coupling amplitude [Pa]: "
                  << coupling.GetAcousticAmplitudePa() << '\n';
        std::cout << "Coupling phase [rad]: "
                  << coupling.GetPhaseRad() << '\n';
        std::cout << "Coupling omega [rad/s]: "
                  << coupling.GetAngularFrequencyRadS() << '\n';
        std::cout << "Coupling static pressure [Pa]: "
                  << coupling.GetStaticPressurePa() << '\n';
        std::cout << "p_inf(0) [Pa]: " << coupling.Evaluate(0.0) << '\n';
        std::cout << "p_inf(T/4) [Pa]: "
                  << coupling.Evaluate(0.25 * period_s) << '\n';
        std::cout << "p_inf(T/2) [Pa]: "
                  << coupling.Evaluate(0.5 * period_s) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error while solving acoustics: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
