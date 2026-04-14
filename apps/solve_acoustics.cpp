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
        const PressureCoupling coupling = PressureCoupling::FromComplexPressure(
            config.liquid.static_pressure_pa,
            problem.Omega,
            pb);
        const double period_s = 2.0 * 3.14159265358979323846 /
                                coupling.GetAngularFrequencyRadS();

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
        std::cout << "Pb phase [rad]: " << std::arg(pb) << '\n';
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
