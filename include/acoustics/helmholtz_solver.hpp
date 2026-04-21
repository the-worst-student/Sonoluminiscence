#pragma once

#include <complex>
#include <string>
#include <vector>

#include "acoustics/acoustics_problem.hpp"
#include "acoustics/boundary_conditions.hpp"
#include "mesh/mesh_tags.hpp"

struct MeshNode {
  double R = 0.0;
  double Z = 0.0;
};

struct TriangleElement {
  int Node0 = -1;
  int Node1 = -1;
  int Node2 = -1;
  int PhysicalTag = -1;
};

struct BoundaryEdge {
  int Node0 = -1;
  int Node1 = -1;
  int PhysicalTag = -1;
};

struct HelmholtzSolution {
  std::vector<std::complex<double>> NodalValues;
};

struct TriangleGeometry {
  MeshNode Node0;
  MeshNode Node1;
  MeshNode Node2;

  double JacobianDeterminant = 0.0;
  double Area = 0.0;
};

struct Gradient2D {
  double Dr = 0.0;
  double Dz = 0.0;
};

struct LocalElementMatrix {
  std::complex<double> Values[3][3] = {};
};

struct LocalEdgeVector {
  std::complex<double> Values[2] = {};
};

class HelmholtzSolver {
 public:
  HelmholtzSolver(
      const AcousticsProblem& problem,
      const BoundaryConditions& boundary_conditions);

  void LoadMesh(const std::string& mesh_path);

  void Assemble();

  void Solve();

  std::complex<double> SamplePressure(double r, double z) const;

  const HelmholtzSolution& GetSolution() const;
  const std::vector<MeshNode>& GetNodes() const;
  const std::vector<TriangleElement>& GetElements() const;

  bool IsMeshLoaded() const;
  bool IsAssembled() const;
  bool IsSolved() const;

 private:
  const AcousticsProblem& Problem;
  const BoundaryConditions& Conditions;

  std::vector<MeshNode> Nodes;
  std::vector<TriangleElement> Elements;
  std::vector<BoundaryEdge> BoundaryEdges;

  std::vector<std::vector<std::complex<double>>> Matrix;
  std::vector<std::complex<double>> RightHandSide;
  HelmholtzSolution Solution;

  bool MeshLoaded = false;
  bool Assembled = false;
  bool Solved = false;

  std::size_t GetNumberOfDofs() const;
  void ValidateMeshData() const;
  void InitializeLinearSystem();
  void ClearSolution();

  TriangleGeometry BuildTriangleGeometry(const TriangleElement& element) const;

  std::vector<Gradient2D> ComputeLinearShapeGradients(
      const TriangleGeometry& geometry) const;

  double ComputeElementCentroidR(const TriangleGeometry& geometry) const;

  LocalElementMatrix ComputeLocalStiffnessMatrix(
      const TriangleGeometry& geometry,
      const std::vector<Gradient2D>& gradients) const;

  LocalElementMatrix ComputeLocalMassMatrix(
      const TriangleGeometry& geometry) const;

  LocalElementMatrix ComputeLocalHelmholtzMatrix(
      const TriangleGeometry& geometry,
      const std::vector<Gradient2D>& gradients) const;

  double ComputeBoundaryEdgeLength(const BoundaryEdge& edge) const;

  double ComputeBoundaryEdgeCentroidR(const BoundaryEdge& edge) const;

  LocalEdgeVector ComputeLocalSourceEdgeVector(
      const BoundaryEdge& edge) const;

  void AddLocalElementMatrixToGlobalSystem(
      const TriangleElement& element,
      const LocalElementMatrix& local_matrix);

  void AddLocalEdgeVectorToGlobalRightHandSide(
      const BoundaryEdge& edge,
      const LocalEdgeVector& local_vector);

  bool IsFluidElement(const TriangleElement& element) const;

  void ValidateAssembledSystem() const;
};
