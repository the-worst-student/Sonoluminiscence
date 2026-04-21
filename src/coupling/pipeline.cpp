#include <coupling/pipeline.hpp>

#include <complex>
#include <sstream>
#include <vector>

#include <acoustics/acoustics_problem.hpp>
#include <acoustics/boundary_conditions.hpp>
#include <acoustics/field_sampler.hpp>
#include <acoustics/helmholtz_solver.hpp>
#include <bubble/ode_solver.hpp>
#include <coupling/pressure_coupling.hpp>
#include <io/bubble_vtk_writer.hpp>
#include <io/csv_writer.hpp>
#include <io/result_writer.hpp>
#include <io/vtk_writer.hpp>
#include <mesh/gmsh_driver.hpp>
#include <postprocess/metrics.hpp>

Pipeline::Pipeline(const ProjectConfig& config) : Config(config) {
}

PipelineResult Pipeline::Run(const PipelineOutputs& outputs) const {
    GmshDriver Driver(Config);
    Driver.BuildAxisymmetricMesh(outputs.MeshPath);

    const AcousticsProblem Problem = AcousticsProblem::FromConfig(Config);
    const BoundaryConditions Conditions(Problem);

    HelmholtzSolver AcousticSolver(Problem, Conditions);
    AcousticSolver.LoadMesh(outputs.MeshPath);
    AcousticSolver.Assemble();
    AcousticSolver.Solve();

    const double BubbleR = Config.geometry.bubble_position.r_m;
    const double BubbleZ = Config.geometry.bubble_position.z_m;

    const std::complex<double> Pb =
        AcousticSolver.SamplePressure(BubbleR, BubbleZ);

    const PressureCoupling Coupling =
        PressureCoupling::FromComplexPressure(
            Config.liquid.static_pressure_pa,
            Problem.Omega,
            Pb);

    const AcousticFieldData FieldData =
        FieldSampler::BuildFromSolver(AcousticSolver);

    VtkWriter::WriteAxisymmetricSlice(
        outputs.SliceVtkPath,
        FieldData,
        BubbleR,
        BubbleZ);

    VtkWriter::WritePseudoVolume(
        outputs.PseudoVolumeVtkPath,
        FieldData,
        BubbleR,
        BubbleZ,
        outputs.PseudoVolumeSectors);

    OdeSolver BubbleSolver(Config, Coupling);
    BubbleSolver.Solve();

    const OdeSolveResult& SolveResult = BubbleSolver.GetResult();
    const std::vector<BubbleSample>& Samples = SolveResult.Samples;

    if (!Samples.empty()) {
        CsvWriter::WriteBubbleSamples(outputs.BubbleCsvPath, Samples);

        BubbleVtkWriter::Options BubbleOptions;
        BubbleOptions.OutputDirectory = outputs.BubbleFramesDirectory;
        BubbleOptions.AzimuthalSegments = outputs.PseudoVolumeSectors;
        BubbleOptions.PolarSegments = 24;

        BubbleVtkWriter::WriteAnimationFrames(
            Samples,
            BubbleR,
            BubbleZ,
            BubbleOptions);
    }

    BubbleMetrics MetricsValue =
        Metrics::EvaluateBubbleMetrics(SolveResult, Config);

    const LuminescenceCriteria Criteria;
    Metrics::ApplyLuminescenceCriteria(MetricsValue, Criteria);

    PipelineResult Result;
    Result.BubblePressurePa = Pb;
    Result.BubblePressureAbsPa = std::abs(Pb);
    Result.BubblePressurePhaseRad = std::arg(Pb);
    Result.SavedSamples = Samples.size();
    Result.Metrics = MetricsValue;
    Result.Outputs = outputs;

    std::ostringstream Summary;
    Summary << "Pipeline finished\n";
    Summary << "Mesh: " << outputs.MeshPath << '\n';
    Summary << "Bubble CSV: " << outputs.BubbleCsvPath << '\n';
    Summary << "Acoustic slice VTK: " << outputs.SliceVtkPath << '\n';
    Summary << "Acoustic pseudovolume VTK: " << outputs.PseudoVolumeVtkPath << '\n';
    Summary << "Bubble animation frames: " << outputs.BubbleFramesDirectory << '\n';

    Summary << "Pb real [Pa]: " << std::real(Pb) << '\n';
    Summary << "Pb imag [Pa]: " << std::imag(Pb) << '\n';
    Summary << "Pb abs [Pa]: " << std::abs(Pb) << '\n';
    Summary << "Pb phase [rad]: " << std::arg(Pb) << '\n';

    Summary << "Coupling amplitude [Pa]: "
            << Coupling.GetAcousticAmplitudePa() << '\n';
    Summary << "Coupling phase [rad]: "
            << Coupling.GetPhaseRad() << '\n';

    Summary << "Bubble solve success: "
            << std::boolalpha << SolveResult.Success << '\n';
    Summary << "Bubble stop reason: "
            << OdeSolver::StopReasonToString(SolveResult.StopReason) << '\n';
    Summary << "Saved bubble samples: " << Result.SavedSamples << '\n';

    Summary << "Initial radius [m]: " << MetricsValue.InitialRadiusM << '\n';
    Summary << "Initial temperature [K]: " << MetricsValue.InitialTemperatureK << '\n';
    Summary << "Static pressure [Pa]: " << MetricsValue.StaticPressurePa << '\n';
    Summary << "Hard-core radius [m]: " << MetricsValue.HardCoreRadiusM << '\n';

    Summary << "Min radius [m]: " << MetricsValue.MinRadiusM << '\n';
    Summary << "Time at min radius [s]: " << MetricsValue.TimeAtMinRadiusS << '\n';

    Summary << "Max gas temperature [K]: " << MetricsValue.MaxTemperatureK << '\n';
    Summary << "Time at max gas temperature [s]: "
            << MetricsValue.TimeAtMaxTemperatureS << '\n';

    Summary << "Max gas pressure [Pa]: " << MetricsValue.MaxGasPressurePa << '\n';
    Summary << "Time at max gas pressure [s]: "
            << MetricsValue.TimeAtMaxGasPressureS << '\n';

    Summary << "Max inward speed [m/s]: "
            << MetricsValue.MaxInwardSpeedMPerS << '\n';
    Summary << "Time at max inward speed [s]: "
            << MetricsValue.TimeAtMaxInwardSpeedS << '\n';

    Summary << "Final time [s]: " << MetricsValue.FinalTimeS << '\n';
    Summary << "Final radius [m]: " << MetricsValue.FinalRadiusM << '\n';
    Summary << "Final temperature [K]: " << MetricsValue.FinalTemperatureK << '\n';
    Summary << "Final gas pressure [Pa]: " << MetricsValue.FinalGasPressurePa << '\n';
    Summary << "Final external pressure [Pa]: "
            << MetricsValue.FinalExternalPressurePa << '\n';

    Summary << "Compression ratio R0/Rmin: "
            << MetricsValue.CompressionRatio << '\n';
    Summary << "Hard-core proximity ratio Rmin/Rhc: "
            << MetricsValue.HardCoreProximityRatio << '\n';
    Summary << "Temperature gain Tmax/T0: "
            << MetricsValue.TemperatureGain << '\n';
    Summary << "Pressure gain Pg,max/P0: "
            << MetricsValue.PressureGain << '\n';

    Summary << "Collapse detected: "
            << std::boolalpha << MetricsValue.CollapseDetected << '\n';
    if (MetricsValue.CollapseDetected) {
        Summary << "Collapse time [s]: " << MetricsValue.CollapseTimeS << '\n';
        Summary << "Radius at collapse [m]: "
                << MetricsValue.RadiusAtCollapseM << '\n';
        Summary << "Temperature at collapse [K]: "
                << MetricsValue.TemperatureAtCollapseK << '\n';
        Summary << "Gas pressure at collapse [Pa]: "
                << MetricsValue.GasPressureAtCollapsePa << '\n';
    }

    Summary << "Strong compression: "
            << std::boolalpha << MetricsValue.StrongCompression << '\n';
    Summary << "Strong heating: "
            << std::boolalpha << MetricsValue.StrongHeating << '\n';
    Summary << "Strong gas pressure: "
            << std::boolalpha << MetricsValue.StrongGasPressure << '\n';
    Summary << "Strong inward speed: "
            << std::boolalpha << MetricsValue.StrongInwardSpeed << '\n';
    Summary << "Near hard-core: "
            << std::boolalpha << MetricsValue.NearHardCore << '\n';
    Summary << "Potential luminescence: "
            << std::boolalpha << MetricsValue.PotentialLuminescence << '\n';

    ResultWriter::WriteTextFile(outputs.SummaryPath, Summary.str());

    return Result;
}