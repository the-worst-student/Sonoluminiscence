#include <acoustics/gorkov_candidates.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {

double SafeNormalize(
    const double value,
    const double min_value,
    const double max_value) {
    const double Range = max_value - min_value;
    if (std::abs(Range) < 1e-30) {
        return 0.0;
    }
    return (value - min_value) / Range;
}

std::vector<std::vector<std::size_t>> BuildNodeAdjacency(
    const AcousticFieldData& field_data) {
    std::vector<std::unordered_set<std::size_t>> Sets(field_data.Nodes.size());

    for (const TriangleElement& Element : field_data.Elements) {
        const std::size_t A = static_cast<std::size_t>(Element.Node0);
        const std::size_t B = static_cast<std::size_t>(Element.Node1);
        const std::size_t C = static_cast<std::size_t>(Element.Node2);
        Sets[A].insert(B);
        Sets[A].insert(C);
        Sets[B].insert(A);
        Sets[B].insert(C);
        Sets[C].insert(A);
        Sets[C].insert(B);
    }

    std::vector<std::vector<std::size_t>> Result(Sets.size());
    for (std::size_t i = 0; i < Sets.size(); ++i) {
        Result[i].assign(Sets[i].begin(), Sets[i].end());
    }
    return Result;
}

bool IsLocalPotentialMinimum(
    const AcousticFieldData& field_data,
    const std::vector<std::vector<std::size_t>>& adjacency,
    const std::size_t node_index) {
    const double Value = field_data.NodalFields[node_index].GorkovPotentialJ;
    for (const std::size_t Neighbor : adjacency[node_index]) {
        if (field_data.NodalFields[Neighbor].GorkovPotentialJ < Value) {
            return false;
        }
    }
    return !adjacency[node_index].empty();
}

std::string CandidateIdString(const int candidate_id) {
    std::ostringstream Output;
    Output << "candidate_" << std::setw(3) << std::setfill('0') << candidate_id;
    return Output.str();
}

BubbleExcitation BuildExcitation(
    const GorkovCandidate& candidate,
    const ProjectConfig& config,
    const double angular_frequency_rad_s) {
    BubbleExcitation Excitation;
    Excitation.CaseId = CandidateIdString(candidate.CandidateId);
    if (!config.project.run_id.empty()) {
        Excitation.CaseId = config.project.run_id + "_" + Excitation.CaseId;
    }
    Excitation.FrequencyHz = config.acoustics.frequency_hz;
    Excitation.AngularFrequencyRadS = angular_frequency_rad_s;
    Excitation.BubbleRM = candidate.RM;
    Excitation.BubbleZM = candidate.ZM;
    Excitation.StaticPressurePa = config.liquid.static_pressure_pa;
    Excitation.LiquidTemperatureK = config.liquid.temperature_k;
    Excitation.EquilibriumRadiusM = config.bubble.equilibrium_radius_m;
    Excitation.BubblePressurePa = {candidate.PressureRealPa, candidate.PressureImagPa};
    return Excitation;
}

}  // namespace

std::vector<GorkovCandidate> GorkovCandidates::SelectCandidates(
    const AcousticFieldData& field_data,
    const GorkovCandidateOptions& options) {
    if (field_data.Nodes.size() != field_data.NodalFields.size()) {
        throw std::runtime_error("Node count does not match field count for candidates");
    }

    if (field_data.Nodes.empty()) {
        return {};
    }

    double MaxPressure = 0.0;
    double MinPotential = std::numeric_limits<double>::max();
    double MaxPotential = -std::numeric_limits<double>::max();
    double MaxForce = 0.0;

    for (const AcousticNodalField& Field : field_data.NodalFields) {
        MaxPressure = std::max(MaxPressure, Field.PressureAbsPa);
        MinPotential = std::min(MinPotential, Field.GorkovPotentialJ);
        MaxPotential = std::max(MaxPotential, Field.GorkovPotentialJ);
        MaxForce = std::max(MaxForce, Field.GorkovForceAbsN);
    }

    const std::vector<std::vector<std::size_t>> Adjacency = BuildNodeAdjacency(field_data);
    std::vector<GorkovCandidate> Pool;

    const double MinPressure = options.MinPressureFraction * MaxPressure;
    for (std::size_t i = 0; i < field_data.Nodes.size(); ++i) {
        const AcousticNodalField& Field = field_data.NodalFields[i];
        if (Field.PressureAbsPa < MinPressure) {
            continue;
        }

        const bool LocalMinimum = IsLocalPotentialMinimum(field_data, Adjacency, i);
        const double PressureNorm = MaxPressure > 0.0 ? Field.PressureAbsPa / MaxPressure : 0.0;
        const double PotentialNorm = SafeNormalize(Field.GorkovPotentialJ, MinPotential, MaxPotential);
        const double ForceNorm = MaxForce > 0.0 ? Field.GorkovForceAbsN / MaxForce : 0.0;
        const double LocalBonus = LocalMinimum ? 0.15 : 0.0;

        GorkovCandidate Candidate;
        Candidate.NodeIndex = i;
        Candidate.RM = field_data.Nodes[i].R;
        Candidate.ZM = field_data.Nodes[i].Z;
        Candidate.PressureRealPa = Field.PressureRealPa;
        Candidate.PressureImagPa = Field.PressureImagPa;
        Candidate.PressureAbsPa = Field.PressureAbsPa;
        Candidate.PressurePhaseRad = Field.PressurePhaseRad;
        Candidate.DriveAmplitudePa = Field.PressureAbsPa;
        Candidate.DrivePhaseRad = -Field.PressurePhaseRad;
        Candidate.GorkovPotentialJ = Field.GorkovPotentialJ;
        Candidate.GorkovForceAbsN = Field.GorkovForceAbsN;
        Candidate.VelocityAbsMS = Field.VelocityAbsMS;
        Candidate.LocalPotentialMinimum = LocalMinimum;
        Candidate.Score =
            options.PressureWeight * PressureNorm -
            options.PotentialWeight * PotentialNorm -
            options.ForceWeight * ForceNorm + LocalBonus;
        Pool.push_back(Candidate);
    }

    std::sort(Pool.begin(), Pool.end(), [](const GorkovCandidate& left,
                                           const GorkovCandidate& right) {
        return left.Score > right.Score;
    });

    if (options.MaxCandidates > 0 &&
        Pool.size() > static_cast<std::size_t>(options.MaxCandidates)) {
        Pool.resize(static_cast<std::size_t>(options.MaxCandidates));
    }

    for (std::size_t i = 0; i < Pool.size(); ++i) {
        Pool[i].CandidateId = static_cast<int>(i) + 1;
    }

    return Pool;
}

void GorkovCandidates::MarkCandidates(
    AcousticFieldData* field_data,
    const std::vector<GorkovCandidate>& candidates) {
    if (field_data == nullptr) {
        throw std::runtime_error("Cannot mark candidates in null field data");
    }

    for (AcousticNodalField& Field : field_data->NodalFields) {
        Field.CandidateMarker = 0.0;
        Field.CandidateScore = 0.0;
    }

    for (const GorkovCandidate& Candidate : candidates) {
        if (Candidate.NodeIndex >= field_data->NodalFields.size()) {
            continue;
        }
        AcousticNodalField& Field = field_data->NodalFields[Candidate.NodeIndex];
        Field.CandidateMarker = static_cast<double>(Candidate.CandidateId);
        Field.CandidateScore = Candidate.Score;
    }

    for (std::size_t i = 0; i < field_data->Elements.size(); ++i) {
        const TriangleElement& Element = field_data->Elements[i];
        const AcousticNodalField& Field0 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node0));
        const AcousticNodalField& Field1 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node1));
        const AcousticNodalField& Field2 =
            field_data->NodalFields.at(static_cast<std::size_t>(Element.Node2));
        field_data->CellFields[i].CandidateScore =
            (Field0.CandidateScore + Field1.CandidateScore + Field2.CandidateScore) / 3.0;
    }
}

std::string GorkovCandidates::FieldCsv(const AcousticFieldData& field_data) {
    std::ostringstream Output;
    Output << std::setprecision(17);
    Output << "node_index,r_m,z_m,pressure_re_pa,pressure_im_pa,pressure_abs_pa,"
           << "pressure_phase_rad,grad_p_abs_pa_per_m,velocity_abs_m_s,"
           << "gorkov_potential_j,gorkov_force_r_n,gorkov_force_z_n,"
           << "gorkov_force_abs_n,candidate_marker,candidate_score\n";

    for (std::size_t i = 0; i < field_data.Nodes.size(); ++i) {
        const MeshNode& Node = field_data.Nodes[i];
        const AcousticNodalField& Field = field_data.NodalFields[i];
        Output << i << ','
               << Node.R << ','
               << Node.Z << ','
               << Field.PressureRealPa << ','
               << Field.PressureImagPa << ','
               << Field.PressureAbsPa << ','
               << Field.PressurePhaseRad << ','
               << Field.GradPressureAbsPaPerM << ','
               << Field.VelocityAbsMS << ','
               << Field.GorkovPotentialJ << ','
               << Field.GorkovForceRN << ','
               << Field.GorkovForceZN << ','
               << Field.GorkovForceAbsN << ','
               << Field.CandidateMarker << ','
               << Field.CandidateScore << '\n';
    }

    return Output.str();
}

std::string GorkovCandidates::CandidatesCsv(
    const std::vector<GorkovCandidate>& candidates) {
    std::ostringstream Output;
    Output << std::setprecision(17);
    Output << "candidate_id,node_index,r_m,z_m,pressure_re_pa,pressure_im_pa,"
           << "pressure_abs_pa,pressure_phase_rad,drive_amplitude_pa,"
           << "drive_phase_rad,gorkov_potential_j,gorkov_force_abs_n,"
           << "velocity_abs_m_s,local_potential_minimum,score\n";

    for (const GorkovCandidate& Candidate : candidates) {
        Output << CandidateIdString(Candidate.CandidateId) << ','
               << Candidate.NodeIndex << ','
               << Candidate.RM << ','
               << Candidate.ZM << ','
               << Candidate.PressureRealPa << ','
               << Candidate.PressureImagPa << ','
               << Candidate.PressureAbsPa << ','
               << Candidate.PressurePhaseRad << ','
               << Candidate.DriveAmplitudePa << ','
               << Candidate.DrivePhaseRad << ','
               << Candidate.GorkovPotentialJ << ','
               << Candidate.GorkovForceAbsN << ','
               << Candidate.VelocityAbsMS << ','
               << (Candidate.LocalPotentialMinimum ? 1 : 0) << ','
               << Candidate.Score << '\n';
    }

    return Output.str();
}

std::string GorkovCandidates::BubbleExcitationsCsv(
    const std::vector<GorkovCandidate>& candidates,
    const ProjectConfig& config,
    const double angular_frequency_rad_s) {
    std::string Result = BubbleExcitation::CsvHeader();
    for (const GorkovCandidate& Candidate : candidates) {
        Result += BuildExcitation(Candidate, config, angular_frequency_rad_s).CsvRow();
    }
    return Result;
}
