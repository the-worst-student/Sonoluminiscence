#include <algorithm>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
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

std::string FrequencyTag(double frequency_hz) {
    std::ostringstream Output;
    Output << "f_" << std::fixed << std::setprecision(0) << frequency_hz << "Hz";
    return Output.str();
}

std::string CsvEscape(const std::string& value) {
    bool NeedsQuotes = false;
    for (const char Symbol : value) {
        if (Symbol == ',' || Symbol == '"' || Symbol == '\n' || Symbol == '\r') {
            NeedsQuotes = true;
            break;
        }
    }

    if (!NeedsQuotes) {
        return value;
    }

    std::string Result;
    Result.reserve(value.size() + 2);
    Result.push_back('"');
    for (const char Symbol : value) {
        if (Symbol == '"') {
            Result += "\"\"";
        } else {
            Result.push_back(Symbol);
        }
    }
    Result.push_back('"');
    return Result;
}

std::vector<double> BuildFrequencyGrid(
    const double min_frequency_hz,
    const double max_frequency_hz,
    const int steps) {
    if (min_frequency_hz <= 0.0 || max_frequency_hz <= 0.0) {
        throw std::runtime_error("Frequency bounds must be positive");
    }

    if (steps <= 0) {
        throw std::runtime_error("Frequency step count must be positive");
    }

    if (steps == 1) {
        return {min_frequency_hz};
    }

    std::vector<double> Frequencies;
    Frequencies.reserve(static_cast<std::size_t>(steps));

    for (int Index = 0; Index < steps; ++Index) {
        const double Alpha = static_cast<double>(Index) / static_cast<double>(steps - 1);
        Frequencies.push_back(min_frequency_hz + Alpha * (max_frequency_hz - min_frequency_hz));
    }

    return Frequencies;
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

struct FrequencyScanRow {
    double FrequencyHz = 0.0;
    double OmegaRadS = 0.0;
    double ManualPbAbsPa = 0.0;
    double ManualDrivePhaseRad = 0.0;
    double MaxPressureAbsPa = 0.0;
    int CandidateCount = 0;
    std::string BestCandidateId;
    std::size_t BestNodeIndex = 0;
    double BestRM = 0.0;
    double BestZM = 0.0;
    double BestPressureAbsPa = 0.0;
    double BestDrivePhaseRad = 0.0;
    double BestGorkovPotentialJ = 0.0;
    double BestGorkovForceAbsN = 0.0;
    double BestCurvatureTraceJPerM2 = 0.0;
    double BestCurvatureLambdaMinJPerM2 = 0.0;
    double BestCurvatureLambdaMaxJPerM2 = 0.0;
    int BestCurvaturePositiveDefinite = 0;
    double BestScore = 0.0;
    std::string FieldCsvPath;
    std::string CandidatesCsvPath;
    std::string SliceVtkPath;
    std::string PseudovolumeVtkPath;
};

void WriteFrequencyScanSummary(
    const std::string& output_path,
    const std::vector<FrequencyScanRow>& rows) {
    std::ofstream Output(output_path);
    if (!Output.is_open()) {
        throw std::runtime_error("Failed to open frequency scan summary: " + output_path);
    }

    Output << std::setprecision(17);
    Output << "frequency_hz,omega_rad_s,manual_pb_abs_pa,manual_drive_phase_rad,"
           << "max_pressure_abs_pa,candidate_count,best_candidate_id,best_node_index,"
           << "best_r_m,best_z_m,best_pressure_abs_pa,best_drive_phase_rad,"
           << "best_gorkov_potential_j,best_gorkov_force_abs_n,"
           << "best_gorkov_curvature_trace_j_m2,"
           << "best_gorkov_curvature_lambda_min_j_m2,"
           << "best_gorkov_curvature_lambda_max_j_m2,"
           << "best_gorkov_curvature_positive_definite,best_score,"
           << "field_csv_path,candidates_csv_path,slice_vtk_path,pseudovolume_vtk_path\n";

    for (const FrequencyScanRow& Row : rows) {
        Output << Row.FrequencyHz << ','
               << Row.OmegaRadS << ','
               << Row.ManualPbAbsPa << ','
               << Row.ManualDrivePhaseRad << ','
               << Row.MaxPressureAbsPa << ','
               << Row.CandidateCount << ','
               << CsvEscape(Row.BestCandidateId) << ','
               << Row.BestNodeIndex << ','
               << Row.BestRM << ','
               << Row.BestZM << ','
               << Row.BestPressureAbsPa << ','
               << Row.BestDrivePhaseRad << ','
               << Row.BestGorkovPotentialJ << ','
               << Row.BestGorkovForceAbsN << ','
               << Row.BestCurvatureTraceJPerM2 << ','
               << Row.BestCurvatureLambdaMinJPerM2 << ','
               << Row.BestCurvatureLambdaMaxJPerM2 << ','
               << Row.BestCurvaturePositiveDefinite << ','
               << Row.BestScore << ','
               << CsvEscape(Row.FieldCsvPath) << ','
               << CsvEscape(Row.CandidatesCsvPath) << ','
               << CsvEscape(Row.SliceVtkPath) << ','
               << CsvEscape(Row.PseudovolumeVtkPath) << '\n';
    }
}

void PrintUsage() {
    std::cerr
        << "Usage: ./scan_frequency <config.yaml> <mesh.msh> <output_dir> "
        << "<f_min_hz> <f_max_hz> <steps> [max_candidates] [pseudovolume_sectors]\n"
        << "Example: ./scan_frequency configs/base.yaml mesh.msh results/frequency_scan 20000 200000 25 10 0\n";
}

}  // namespace

int main(const int argc, char** argv) {
    try {
        if (argc < 7) {
            PrintUsage();
            return 1;
        }

        const std::string ConfigPath = ResolveExistingPath(argv[1]);
        const std::string MeshPath = ResolveExistingPath(argv[2]);
        const std::filesystem::path OutputDir(argv[3]);
        const double MinFrequencyHz = std::stod(argv[4]);
        const double MaxFrequencyHz = std::stod(argv[5]);
        const int Steps = std::stoi(argv[6]);
        const int MaxCandidates = argc > 7 ? std::stoi(argv[7]) : 10;
        const int PseudovolumeSectors = argc > 8 ? std::stoi(argv[8]) : 0;

        std::filesystem::create_directories(OutputDir);

        const ProjectConfig BaseConfig = YamlReader::ReadProjectConfig(ConfigPath);
        const std::vector<double> Frequencies = BuildFrequencyGrid(
            MinFrequencyHz,
            MaxFrequencyHz,
            Steps);

        std::vector<FrequencyScanRow> Rows;
        Rows.reserve(Frequencies.size());

        std::string AllExcitationsCsv = BubbleExcitation::CsvHeader();

        std::cout << "Frequency scan started\n";
        std::cout << "Config: " << ConfigPath << '\n';
        std::cout << "Mesh: " << MeshPath << '\n';
        std::cout << "Output dir: " << OutputDir.string() << '\n';
        std::cout << "Frequencies: " << Frequencies.size() << '\n';

        for (const double FrequencyHz : Frequencies) {
            ProjectConfig Config = BaseConfig;
            Config.acoustics.frequency_hz = FrequencyHz;
            Config.project.run_id = BaseConfig.project.run_id + "_" + FrequencyTag(FrequencyHz);

            const std::filesystem::path CaseDir = OutputDir / FrequencyTag(FrequencyHz);
            std::filesystem::create_directories(CaseDir);

            const AcousticsProblem Problem = AcousticsProblem::FromConfig(Config);
            const BoundaryConditions Conditions(Problem);
            HelmholtzSolver Solver(Problem, Conditions);

            std::cout << "frequency_hz=" << FrequencyHz << " load/assemble/solve\n";
            Solver.LoadMesh(MeshPath);
            Solver.Assemble();
            Solver.Solve();

            const double BubbleR = Config.geometry.bubble_position.r_m;
            const double BubbleZ = Config.geometry.bubble_position.z_m;
            const std::complex<double> ManualPb = Solver.SamplePressure(BubbleR, BubbleZ);
            const BubbleExcitation ManualExcitation = BuildManualBubbleExcitation(
                Config,
                Problem,
                BubbleR,
                BubbleZ,
                ManualPb);

            AcousticFieldData FieldData = FieldSampler::BuildFromSolver(Solver);
            FieldDerivatives::AddPressureGradient(&FieldData);
            const GorkovParameters GorkovParameters =
                GorkovPotential::FromConfig(Config, Problem.Omega);
            GorkovPotential::AddGorkovFields(&FieldData, GorkovParameters);

            GorkovCandidateOptions CandidateOptions;
            CandidateOptions.MaxCandidates = MaxCandidates;
            const std::vector<GorkovCandidate> Candidates =
                GorkovCandidates::SelectCandidates(FieldData, CandidateOptions);
            GorkovCandidates::MarkCandidates(&FieldData, Candidates);

            const std::filesystem::path FieldCsvPath = CaseDir / "gorkov_field.csv";
            const std::filesystem::path CandidatesCsvPath = CaseDir / "gorkov_candidates.csv";
            const std::filesystem::path SliceVtkPath = CaseDir / "acoustic_slice.vtk";
            const std::filesystem::path PseudovolumeVtkPath = CaseDir / "acoustic_pseudovolume.vtk";

            ResultWriter::WriteTextFile(
                FieldCsvPath.string(),
                GorkovCandidates::FieldCsv(FieldData));
            ResultWriter::WriteTextFile(
                CandidatesCsvPath.string(),
                GorkovCandidates::CandidatesCsv(Candidates));
            VtkWriter::WriteAxisymmetricSlice(
                SliceVtkPath.string(),
                FieldData,
                BubbleR,
                BubbleZ);

            std::string PseudovolumePathString;
            if (PseudovolumeSectors >= 3) {
                VtkWriter::WritePseudoVolume(
                    PseudovolumeVtkPath.string(),
                    FieldData,
                    BubbleR,
                    BubbleZ,
                    PseudovolumeSectors);
                PseudovolumePathString = PseudovolumeVtkPath.string();
            }

            AllExcitationsCsv += ManualExcitation.CsvRow();
            const std::string CandidateExcitationsCsv =
                GorkovCandidates::BubbleExcitationsCsv(Candidates, Config, Problem.Omega);
            AllExcitationsCsv += CandidateExcitationsCsv.substr(BubbleExcitation::CsvHeader().size());

            FrequencyScanRow Row;
            Row.FrequencyHz = FrequencyHz;
            Row.OmegaRadS = Problem.Omega;
            Row.ManualPbAbsPa = std::abs(ManualPb);
            Row.ManualDrivePhaseRad = ManualExcitation.PressurePhaseRad();
            Row.CandidateCount = static_cast<int>(Candidates.size());
            Row.FieldCsvPath = FieldCsvPath.string();
            Row.CandidatesCsvPath = CandidatesCsvPath.string();
            Row.SliceVtkPath = SliceVtkPath.string();
            Row.PseudovolumeVtkPath = PseudovolumePathString;

            for (const AcousticNodalField& Field : FieldData.NodalFields) {
                Row.MaxPressureAbsPa = std::max(Row.MaxPressureAbsPa, Field.PressureAbsPa);
            }

            if (!Candidates.empty()) {
                const GorkovCandidate& Best = Candidates.front();
                Row.BestCandidateId = Config.project.run_id + "_candidate_001";
                Row.BestNodeIndex = Best.NodeIndex;
                Row.BestRM = Best.RM;
                Row.BestZM = Best.ZM;
                Row.BestPressureAbsPa = Best.PressureAbsPa;
                Row.BestDrivePhaseRad = Best.DrivePhaseRad;
                Row.BestGorkovPotentialJ = Best.GorkovPotentialJ;
                Row.BestGorkovForceAbsN = Best.GorkovForceAbsN;
                Row.BestCurvatureTraceJPerM2 = Best.GorkovCurvatureTraceJPerM2;
                Row.BestCurvatureLambdaMinJPerM2 = Best.GorkovCurvatureLambdaMinJPerM2;
                Row.BestCurvatureLambdaMaxJPerM2 = Best.GorkovCurvatureLambdaMaxJPerM2;
                Row.BestCurvaturePositiveDefinite = Best.GorkovCurvaturePositiveDefinite ? 1 : 0;
                Row.BestScore = Best.Score;
            }

            Rows.push_back(Row);

            std::cout << "frequency_hz=" << FrequencyHz
                      << " manual_pb_abs_pa=" << Row.ManualPbAbsPa
                      << " best_pressure_abs_pa=" << Row.BestPressureAbsPa
                      << " best_curvature_lambda_min=" << Row.BestCurvatureLambdaMinJPerM2
                      << " candidates=" << Row.CandidateCount
                      << '\n';
        }

        const std::filesystem::path SummaryPath = OutputDir / "frequency_scan_summary.csv";
        const std::filesystem::path ExcitationsPath = OutputDir / "bubble_excitations_all.csv";
        WriteFrequencyScanSummary(SummaryPath.string(), Rows);
        ResultWriter::WriteTextFile(ExcitationsPath.string(), AllExcitationsCsv);

        std::cout << "Frequency scan finished\n";
        std::cout << "Summary CSV: " << SummaryPath.string() << '\n';
        std::cout << "Bubble excitations CSV: " << ExcitationsPath.string() << '\n';

        return 0;
    } catch (const std::exception& Error) {
        std::cerr << "Error while scanning frequencies: " << Error.what() << '\n';
        return 1;
    }
}
