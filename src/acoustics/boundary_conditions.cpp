#include <acoustics/boundary_conditions.hpp>
#include <mesh/mesh_tags.hpp>

BoundaryConditions::BoundaryConditions(const AcousticsProblem& problem) : Problem(problem) {}
 bool BoundaryConditionInfo::IsKnown() const {
     return Type != BoundaryConditionType::Unknown;
 }
 bool BoundaryConditionInfo::IsRigid() const {
     return Type == BoundaryConditionType::Rigid;
 }
 bool BoundaryConditionInfo::IsSource() const {
     return Type == BoundaryConditionType::Source;
 }
 bool BoundaryConditionInfo::IsAxis() const {
     return Type == BoundaryConditionType::Axis;
 }

bool BoundaryConditions::IsRigid(const int tag) const {
    return Problem.IsRigidBoundary(tag);
}

bool BoundaryConditions::IsSource(const int tag) const {
    return tag == Problem.SourceBoundaryTag;
}

bool BoundaryConditions::IsAxis(const int tag)  {
    return tag == static_cast<int>(BoundaryTag::cAxis);
}

BoundaryConditionInfo BoundaryConditions::DescribeBoundary(int tag) const {
    BoundaryConditionInfo info;

    if (IsSource(tag)) {
        info.Type = BoundaryConditionType::Source;
        info.Flux = Problem.SourceFlux;
        return info;
    }

    if (IsAxis(tag)) {
        info.Type = BoundaryConditionType::Axis;
        info.Flux = std::complex<double>(0.0, 0.0);
        return info;
    }

    if (IsRigid(tag)) {
        info.Type = BoundaryConditionType::Rigid;
        info.Flux = std::complex<double>(0.0, 0.0);
        return info;
    }

    return info;
}