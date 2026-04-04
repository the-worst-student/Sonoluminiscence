#pragma once

#include <complex>

#include "acoustics/acoustics_problem.hpp"

enum class BoundaryConditionType {
    Unknown,
    Rigid,
    Source,
    Axis
};

struct BoundaryConditionInfo {
    BoundaryConditionType Type = BoundaryConditionType::Unknown;
    std::complex<double> Flux = {0, 0};

    bool IsKnown() const;
    bool IsRigid() const;
    bool IsSource() const;
    bool IsAxis() const;
};

class BoundaryConditions {
public:
    explicit BoundaryConditions(const AcousticsProblem& problem);

    BoundaryConditionInfo DescribeBoundary(int tag) const;

    bool IsRigid(int tag) const;
    bool IsSource(int tag) const;
    static bool IsAxis(int tag) ;

private:
    const AcousticsProblem& Problem;
};