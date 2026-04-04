#include "acoustics/helmholtz_solver.hpp"


#include <cmath>
#include <stdexcept>

HelmholtzSolver::HelmholtzSolver(
    const AcousticsProblem& problem,
    const BoundaryConditions& boundary_conditions)
    : Problem(problem), Conditions(boundary_conditions) {}

void HelmholtzSolver::LoadMesh(const std::string& mesh_path) {
  if (mesh_path.empty()) {
    throw std::invalid_argument("Mesh path must not be empty");
  }

  Nodes.clear();
  Elements.clear();
  BoundaryEdges.clear();

  Matrix.clear();
  RightHandSide.clear();
  ClearSolution();

  MeshLoaded = true;
  Assembled = false;
  Solved = false;
}

void HelmholtzSolver::Assemble() {
  if (!MeshLoaded) {
    throw std::runtime_error("Cannot assemble before mesh is loaded");
  }

  ValidateMeshData();
  InitializeLinearSystem();
  ClearSolution();

  for (const TriangleElement& element : Elements) {
      if (!IsFluidElement(element)) {
          continue;
      }

      const TriangleGeometry& geometry = BuildTriangleGeometry(element);
      const std::vector<Gradient2D> gradients = ComputeLinearShapeGradients(geometry);
      const LocalElementMatrix local_matrix = ComputeLocalHelmholtzMatrix(geometry, gradients);

      AddLocalElementMatrixToGlobalSystem(element, local_matrix);
  }

  for (const BoundaryEdge& edge : BoundaryEdges) {
      if (!Conditions.IsSource(edge.PhysicalTag)) {
          continue;
      }

      const LocalEdgeVector local_vector = ComputeLocalSourceEdgeVector(edge);

      AddLocalEdgeVectorToGlobalRightHandSide(edge, local_vector);

      }

  Assembled = true;
  Solved = false;
}

void HelmholtzSolver::Solve() {
  if (!Assembled) {
    throw std::runtime_error("Cannot solve before system is assembled");
  }

  Solution.NodalValues.assign(
      GetNumberOfDofs(), std::complex<double>(0.0, 0.0));

  Solved = true;
}

std::complex<double> HelmholtzSolver::SamplePressure(
    const double r, const double z) const {
  if (!Solved) {
    throw std::runtime_error("Cannot sample pressure before solve");
  }

  (void)r;
  (void)z;
  return std::complex<double>(0.0, 0.0);
}

const HelmholtzSolution& HelmholtzSolver::GetSolution() const {
  if (!Solved) {
    throw std::runtime_error("Solution is not available before solve");
  }

  return Solution;
}

bool HelmholtzSolver::IsMeshLoaded() const {
  return MeshLoaded;
}

bool HelmholtzSolver::IsAssembled() const {
  return Assembled;
}

bool HelmholtzSolver::IsSolved() const {
  return Solved;
}

std::size_t HelmholtzSolver::GetNumberOfDofs() const {
  return Nodes.size();
}

void HelmholtzSolver::ValidateMeshData() const {
  if (Nodes.empty()) {
    throw std::runtime_error("Mesh does not contain nodes");
  }

  if (Elements.empty()) {
    throw std::runtime_error("Mesh does not contain volume elements");
  }

  const int max_node_index = static_cast<int>(Nodes.size()) - 1;

  for (const TriangleElement& element : Elements) {
    if (element.Node0 < 0 || element.Node0 > max_node_index ||
        element.Node1 < 0 || element.Node1 > max_node_index ||
        element.Node2 < 0 || element.Node2 > max_node_index) {
      throw std::runtime_error(
          "Triangle element contains invalid node index");
    }
  }

  for (const BoundaryEdge& edge : BoundaryEdges) {
    if (edge.Node0 < 0 || edge.Node0 > max_node_index ||
        edge.Node1 < 0 || edge.Node1 > max_node_index) {
      throw std::runtime_error(
          "Boundary edge contains invalid node index");
    }
  }
}

void HelmholtzSolver::InitializeLinearSystem() {
  const std::size_t dof_count = GetNumberOfDofs();

  Matrix.assign(
      dof_count,
      std::vector<std::complex<double>>(
          dof_count, std::complex<double>(0.0, 0.0)));

  RightHandSide.assign(dof_count, std::complex<double>(0.0, 0.0));
}

void HelmholtzSolver::ClearSolution() {
  Solution.NodalValues.clear();
}

TriangleGeometry HelmholtzSolver::BuildTriangleGeometry(
    const TriangleElement& element) const {
  const MeshNode& node0 = Nodes.at(static_cast<std::size_t>(element.Node0));
  const MeshNode& node1 = Nodes.at(static_cast<std::size_t>(element.Node1));
  const MeshNode& node2 = Nodes.at(static_cast<std::size_t>(element.Node2));

  const double j11 = node1.R - node0.R;
  const double j12 = node2.R - node0.R;
  const double j21 = node1.Z - node0.Z;
  const double j22 = node2.Z - node0.Z;

  const double determinant = j11 * j22 - j12 * j21;
  const double area = 0.5 * std::abs(determinant);

  if (area <= 0.0) {
    throw std::runtime_error("Triangle element has non-positive area");
  }

  TriangleGeometry geometry;
  geometry.Node0 = node0;
  geometry.Node1 = node1;
  geometry.Node2 = node2;
  geometry.JacobianDeterminant = determinant;
  geometry.Area = area;

  return geometry;
}

std::vector<Gradient2D> HelmholtzSolver::ComputeLinearShapeGradients(
    const TriangleGeometry& geometry) const {
  const double determinant = geometry.JacobianDeterminant;

  if (std::abs(determinant) < 1e-14) {
    throw std::runtime_error("Triangle Jacobian determinant is zero");
  }

  const double inv_det = 1.0 / determinant;

  const double r0 = geometry.Node0.R;
  const double z0 = geometry.Node0.Z;
  const double r1 = geometry.Node1.R;
  const double z1 = geometry.Node1.Z;
  const double r2 = geometry.Node2.R;
  const double z2 = geometry.Node2.Z;

  std::vector<Gradient2D> gradients(3);

  gradients[0].Dr = (z1 - z2) * inv_det;
  gradients[0].Dz = (r2 - r1) * inv_det;

  gradients[1].Dr = (z2 - z0) * inv_det;
  gradients[1].Dz = (r0 - r2) * inv_det;

  gradients[2].Dr = (z0 - z1) * inv_det;
  gradients[2].Dz = (r1 - r0) * inv_det;

  return gradients;
}

double HelmholtzSolver::ComputeElementCentroidR(
    const TriangleGeometry& geometry) const {
  return (geometry.Node0.R + geometry.Node1.R + geometry.Node2.R) / 3.0;
}

LocalElementMatrix HelmholtzSolver::ComputeLocalStiffnessMatrix(
    const TriangleGeometry& geometry,
    const std::vector<Gradient2D>& gradients) const {
  if (gradients.size() != 3) {
    throw std::runtime_error(
        "Linear triangle must have exactly 3 gradients");
  }

  const double centroid_r = ComputeElementCentroidR(geometry);
  const double scale = centroid_r * geometry.Area;

  LocalElementMatrix matrix;

  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      const double dot_product =
          gradients[col].Dr * gradients[row].Dr +
          gradients[col].Dz * gradients[row].Dz;

      matrix.Values[row][col] =
          std::complex<double>(scale * dot_product, 0.0);
    }
  }

  return matrix;
}

LocalElementMatrix HelmholtzSolver::ComputeLocalMassMatrix(
    const TriangleGeometry& geometry) const {
  const double centroid_r = ComputeElementCentroidR(geometry);
  const double scale = centroid_r * geometry.Area / 12.0;

  LocalElementMatrix matrix;

  matrix.Values[0][0] = std::complex<double>(2.0 * scale, 0.0);
  matrix.Values[0][1] = std::complex<double>(1.0 * scale, 0.0);
  matrix.Values[0][2] = std::complex<double>(1.0 * scale, 0.0);

  matrix.Values[1][0] = std::complex<double>(1.0 * scale, 0.0);
  matrix.Values[1][1] = std::complex<double>(2.0 * scale, 0.0);
  matrix.Values[1][2] = std::complex<double>(1.0 * scale, 0.0);

  matrix.Values[2][0] = std::complex<double>(1.0 * scale, 0.0);
  matrix.Values[2][1] = std::complex<double>(1.0 * scale, 0.0);
  matrix.Values[2][2] = std::complex<double>(2.0 * scale, 0.0);

  return matrix;
}

LocalElementMatrix HelmholtzSolver::ComputeLocalHelmholtzMatrix(
    const TriangleGeometry& geometry,
    const std::vector<Gradient2D>& gradients) const {
  const LocalElementMatrix stiffness =
      ComputeLocalStiffnessMatrix(geometry, gradients);
  const LocalElementMatrix mass =
      ComputeLocalMassMatrix(geometry);

  LocalElementMatrix matrix;

  const double k_squared = Problem.WaveNumber * Problem.WaveNumber;

  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      matrix.Values[row][col] =
          stiffness.Values[row][col] -
          std::complex<double>(k_squared, 0.0) * mass.Values[row][col];
    }
  }

  return matrix;
}

double HelmholtzSolver::ComputeBoundaryEdgeLength(const BoundaryEdge& edge) const {
    const MeshNode& node0 = Nodes.at(static_cast<std::size_t>(edge.Node0));
    const MeshNode& node1 = Nodes.at(static_cast<std::size_t>(edge.Node1));

    const double dr = node1.R - node0.R;
    const double dz = node1.Z - node0.Z;

    return std::sqrt(dr * dr + dz * dz);
}

double HelmholtzSolver::ComputeBoundaryEdgeCentroidR(const BoundaryEdge& edge) const {
    const MeshNode& node0 = Nodes.at(static_cast<std::size_t>(edge.Node0));
    const MeshNode& node1 = Nodes.at(static_cast<std::size_t>(edge.Node1));

    return (node0.R + node1.R) / 2.0;
}

LocalEdgeVector HelmholtzSolver::ComputeLocalSourceEdgeVector(const BoundaryEdge& edge) const {
    const double edge_length = ComputeBoundaryEdgeLength(edge);

    if (edge_length <= 0.0) {
        throw std::runtime_error("Boundary edge length has non-positive length");
    }

    const double centroid_r = ComputeBoundaryEdgeCentroidR(edge);
    const std::complex<double> scale = Problem.SourceFlux * centroid_r * edge_length * 0.5;

    LocalEdgeVector vector;
    vector.Values[0] = scale;
    vector.Values[1] = scale;

    return vector;
}

bool HelmholtzSolver::IsFluidElement(const TriangleElement& element) const{
    return element.PhysicalTag == static_cast<int>(SurfaceTag::cFluid);
}

void HelmholtzSolver::AddLocalElementMatrixToGlobalSystem(const TriangleElement& element, const LocalElementMatrix& local_matrix){
    const int global_indices[3] = {element.Node0, element.Node1, element.Node2};

    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            Matrix[static_cast<std::size_t>(global_indices[row])][static_cast<std::size_t>(col)] += local_matrix.Values[row][col];
        }
    }
}

void HelmholtzSolver::AddLocalEdgeVectorToGlobalRightHandSide(const BoundaryEdge& edge, const LocalEdgeVector& local_vector) {
    const int global_indices[2] = {edge.Node0, edge.Node1};

    for (int row = 0; row < 2; ++row) {
        RightHandSide[static_cast<std::size_t>(global_indices[row])] +=
            local_vector.Values[row];
    }
}

void HelmholtzSolver::ValidateAssembledSystem() const {
    const std::size_t dof_count = GetNumberOfDofs();

    if (Matrix.size() != dof_count) {
        throw std::runtime_error("Matrix size does not match number of DOFs");
    }

    if (RightHandSide.size() != dof_count) {
        throw std::runtime_error("RightHandSide size does not match number of DOFs");
    }

    bool matrix_nonzero = false;
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t col = 0; col < 3; ++col) {
            if (std::abs(Matrix[row][col]) > 1e-14) {
                matrix_nonzero = true;
                break;
            }
        }
        if (matrix_nonzero) {
            break;
        }
    }
}