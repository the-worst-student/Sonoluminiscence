#pragma once

#include <string>
#include <vector>

#include <acoustics/field_sampler.hpp>
#include <coupling/bubble_excitation.hpp>
#include <core/config.hpp>

struct GorkovCandidate {
    int CandidateId = 0;
    std::size_t NodeIndex = 0;
    double RM = 0.0;
    double ZM = 0.0;
    double PressureRealPa = 0.0;
    double PressureImagPa = 0.0;
    double PressureAbsPa = 0.0;
    double PressurePhaseRad = 0.0;
    double DriveAmplitudePa = 0.0;
    double DrivePhaseRad = 0.0;
    double GorkovPotentialJ = 0.0;
    double GorkovForceAbsN = 0.0;
    double VelocityAbsMS = 0.0;
    double Score = 0.0;
    bool LocalPotentialMinimum = false;
};

struct GorkovCandidateOptions {
    int MaxCandidates = 10;
    double MinPressureFraction = 0.05;
    double PressureWeight = 1.0;
    double PotentialWeight = 0.6;
    double ForceWeight = 0.4;
};

class GorkovCandidates {
public:
    static std::vector<GorkovCandidate> SelectCandidates(
        const AcousticFieldData& field_data,
        const GorkovCandidateOptions& options);

    static void MarkCandidates(
        AcousticFieldData* field_data,
        const std::vector<GorkovCandidate>& candidates);

    static std::string FieldCsv(const AcousticFieldData& field_data);

    static std::string CandidatesCsv(
        const std::vector<GorkovCandidate>& candidates);

    static std::string BubbleExcitationsCsv(
        const std::vector<GorkovCandidate>& candidates,
        const ProjectConfig& config,
        double angular_frequency_rad_s);
};
