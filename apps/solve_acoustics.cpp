#include <complex>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "acoustics/acoustics_problem.hpp"
#include "acoustics/boundary_conditions.hpp"
#include "acoustics/field_derivatives.hpp"
#include "acoustics/field_sampler.hpp"
#include "acoustics/gorkov_candidates.hpp"
#include "acoustics/gorkov_potential.hpp"
#include "acoustics/helmholtz_solver.hpp"
#include "coupling/bubble_excitation.hpp"
#include "coupling/pressure_coupling.hpp"
#include "io/result_writer.hpp"
#include "io/vtk_writer.hpp"
#include "io/yaml_reader.hpp"

namespace {

std::string ResolveExistingPath(const std::string& path) {
    namespace fs = std::filesystem;

    const fs::path Direct(path);
    if (fs::exists(Direct)) {
        return Direct.string();
    }

    const fs::path Parent = fs::path("..") / Direct;
    if (fs::exists(Parent)) {
        return Parent.string();
    }

    return path;
}

BubbleExcitation BuildManualBubbleExcitation(
    const ProjectConfig& config,
    const AcousticsProblem& problem,
    const double bubble_r,
    const double bubble_z,
    const std::complex<double>& pb) {
    BubbleExcitation Excitation;
    Excitation.CaseId =
        config.project.run_id.empty() ? "manual_config_point" :
        config.project.run_id + "_manual_config_point";
    Excitation.FrequencyHz = config.acoustics.frequency_hz;
    Excitation.AngularFrequencyRadS = problem.Omega;
    Excitation.BubbleRM = bubble_r;
    Excitation.BubbleZM = bubble_z;
    Excitation.StaticPressurePa = config.liquid.static_pressure_pa;
    Excitation.LiquidTemperatureK = config.liquid.temperature_k;
    Excitation.EquilibriumRadiusM = config.bubble.equilibrium_radius_m;
    Excitation.BubblePressurePa = pb;
    return Excitation;
}

}  // namespace

int main(const int argc, char** argv) {
    const std::string config_path =
        argc > 1 ? ResolveExistingPath(argv[1]) : ResolveExistingPath("configs/base.yaml");
    const std::string mesh_path =
        argc > 2 ? ResolveExistingPath(argv[2]) : ResolveExistingPath("mesh.msh");
    const std::string excitation_path =
        argc > 3 ? argv[3] : "bubble_excitations.csv";
    const std::string slice_vtk_path =
        argc > 4 ? argv[4] : "acoustic_slice.vtk";
    const std::string pseudovolume_vtk_path =
        argc > 5 ? argv[5] : "acoustic_pseudovolume.vtk";
    const int pseudovolume_sectors =
        argc > 6 ? std::stoi(argv[6]) : 48;
    const std::string gorkov_field_path =
        argc > 7 ? argv[7] : "gorkov_field.csv";
    const std::string gorkov_candidates_path =
        argc > 8 ? argv[8] : "gorkov_candidates.csv";
    const int max_candidates =
        argc > 9 ? std::stoi(argv[9]) : 10;

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

        const BubbleExcitation manual_excitation = BuildManualBubbleExcitation(
            config,
            problem,
            bubble_r,
            bubble_z,
            pb);

        const PressureCoupling coupling(
            manual_excitation.StaticPressurePa,
            manual_excitation.PressureAmplitudePa(),
            manual_excitation.AngularFrequencyRadS,
            manual_excitation.PressurePhaseRad());
        const double period_s =
            2.0 * 3.14159265358979323846 / coupling.GetAngularFrequencyRadS();

        AcousticFieldData field_data = FieldSampler::BuildFromSolver(solver);
        FieldDerivatives::AddPressureGradient(&field_data);

        const GorkovParameters gorkov_parameters =
            GorkovPotential::FromConfig(config, problem.Omega);
        GorkovPotential::AddGorkovFields(&field_data, gorkov_parameters);

        GorkovCandidateOptions candidate_options;
        candidate_options.MaxCandidates = max_candidates;
        std::vector<GorkovCandidate> candidates =
            GorkovCandidates::SelectCandidates(field_data, candidate_options);
        GorkovCandidates::MarkCandidates(&field_data, candidates);

        std::string excitation_csv = BubbleExcitation::CsvHeader();
        excitation_csv += manual_excitation.CsvRow();
        excitation_csv += GorkovCandidates::BubbleExcitationsCsv(
            candidates,
            config,
            problem.Omega).substr(BubbleExcitation::CsvHeader().size());
        ResultWriter::WriteTextFile(excitation_path, excitation_csv);

        ResultWriter::WriteTextFile(
            gorkov_field_path,
            GorkovCandidates::FieldCsv(field_data));
        ResultWriter::WriteTextFile(
            gorkov_candidates_path,
            GorkovCandidates::CandidatesCsv(candidates));

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

        std::cout << "Manual bubble point r [m]: " << bubble_r << '\n';
        std::cout << "Manual bubble point z [m]: " << bubble_z << '\n';
        std::cout << "Manual Pb real [Pa]: " << std::real(pb) << '\n';
        std::cout << "Manual Pb imag [Pa]: " << std::imag(pb) << '\n';
        std::cout << "Manual Pb abs [Pa]: " << std::abs(pb) << '\n';
        std::cout << "Manual Pb argument [rad]: "
                  << manual_excitation.PressureArgumentRad() << '\n';
        std::cout << "Manual drive phase [rad]: "
                  << manual_excitation.PressurePhaseRad() << '\n';

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

        std::cout << "Gorkov candidates: " << candidates.size() << '\n';
        if (!candidates.empty()) {
            const GorkovCandidate& Best = candidates.front();
            std::cout << "Best candidate r [m]: " << Best.RM << '\n';
            std::cout << "Best candidate z [m]: " << Best.ZM << '\n';
            std::cout << "Best candidate pressure abs [Pa]: "
                      << Best.PressureAbsPa << '\n';
            std::cout << "Best candidate drive phase [rad]: "
                      << Best.DrivePhaseRad << '\n';
            std::cout << "Best candidate U_G [J]: "
                      << Best.GorkovPotentialJ << '\n';
            std::cout << "Best candidate |F_G| [N]: "
                      << Best.GorkovForceAbsN << '\n';
            std::cout << "Best candidate score: " << Best.Score << '\n';
        }

        std::cout << "Bubble excitations CSV: " << excitation_path << '\n';
        std::cout << "Gorkov field CSV: " << gorkov_field_path << '\n';
        std::cout << "Gorkov candidates CSV: " << gorkov_candidates_path << '\n';
        std::cout << "Acoustic slice VTK: " << slice_vtk_path << '\n';
        std::cout << "Acoustic pseudovolume VTK: "
                  << pseudovolume_vtk_path << '\n';
        std::cout << "Pseudo-volume sectors: "
                  << pseudovolume_sectors << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error while solving acoustics: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
