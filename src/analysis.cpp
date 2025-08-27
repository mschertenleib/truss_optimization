#include "analysis.hpp"

#include <CDT.h>

#include <Eigen/SparseCholesky>

#include <cassert>
#include <cmath>
#include <iostream>

namespace
{

[[nodiscard]] constexpr const char *
to_string(Eigen::ComputationInfo computation_info) noexcept
{
    switch (computation_info)
    {
    case Eigen::Success: return "Success";
    case Eigen::NumericalIssue: return "NumericalIssue";
    case Eigen::NoConvergence: return "NoConvergence";
    case Eigen::InvalidInput: return "InvalidInput";
    }
    return "UNDEFINED";
}

} // namespace

void assemble(Analysis_state &state)
{
    constexpr float young_modulus {200e9f};
    constexpr float area {0.000025f};

    const auto num_dofs = state.nodes.size() * 2;
    const auto num_free_dofs = num_dofs - state.fixed_dofs.size();

    std::vector<int> dof_all_to_free(num_dofs);
    int free_dof {0};
    std::size_t fixed_dof_index {0};
    for (std::uint32_t i {0}; i < num_dofs; ++i)
    {
        if (fixed_dof_index < state.fixed_dofs.size() &&
            i == state.fixed_dofs[fixed_dof_index])
        {
            dof_all_to_free[i] = -1;
            ++fixed_dof_index;
        }
        else
        {
            dof_all_to_free[i] = free_dof;
            ++free_dof;
        }
    }

    std::vector<Eigen::Triplet<float>> triplets;

    state.stiffness_constants.resize(state.elements.size());
    state.element_directions.resize(state.elements.size());

    assert(state.elements.size() == state.activations.size());

    for (std::size_t element_index {0}; element_index < state.elements.size();
         ++element_index)
    {
        const auto [node_i, node_j] = state.elements[element_index];
        const auto activation = state.activations[element_index];

        const auto vec = state.nodes[node_j] - state.nodes[node_i];
        const auto length = norm(vec);
        const auto dir = vec / length;
        state.element_directions[element_index] = dir;
        const auto EA_over_L = activation * (young_modulus * area) / length;
        state.stiffness_constants[element_index] = EA_over_L;

        const auto [c, s] = dir;
        const float element_stiffness_matrix[4][4] {
            {c * c, c * s, -c * c, -c * s},
            {c * s, s * s, -c * s, -s * s},
            {-c * c, -c * s, c * c, c * s},
            {-c * s, -s * s, c * s, s * s}};

        const int dof_indices[4] {dof_all_to_free[2 * node_i],
                                  dof_all_to_free[2 * node_i + 1],
                                  dof_all_to_free[2 * node_j],
                                  dof_all_to_free[2 * node_j + 1]};
        for (unsigned int a {0}; a < 4; ++a)
        {
            for (unsigned int b {0}; b < 4; ++b)
            {
                if (dof_indices[a] != -1 && dof_indices[b] != -1)
                {
                    triplets.emplace_back(dof_indices[a],
                                          dof_indices[b],
                                          EA_over_L *
                                              element_stiffness_matrix[a][b]);
                }
            }
        }
    }

    state.stiffness_matrix.setZero();
    state.stiffness_matrix.resize(static_cast<int>(num_free_dofs),
                                  static_cast<int>(num_free_dofs));
    state.stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    state.stiffness_matrix.prune(0.0f);
}

Linear_system_result solve(const Analysis_state &state)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;
    solver.compute(state.stiffness_matrix);
    solver.factorize(state.stiffness_matrix);

    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Decomposition failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    const auto num_dofs = state.nodes.size() * 2;
    const auto num_free_dofs = num_dofs - state.fixed_dofs.size();
    Eigen::VectorXf free_displacements(num_free_dofs);
    free_displacements = solver.solve(state.loads);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Solving failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    Linear_system_result result {};
    result.displacements.setZero(num_dofs);
    result.displacements(state.free_dofs) = free_displacements;

    result.axial_forces = state.stiffness_constants;
    result.energies = state.stiffness_constants;
    for (std::size_t element_index {0}; element_index < state.elements.size();
         ++element_index)
    {
        const auto [node_i, node_j] = state.elements[element_index];
        const vec2 relative_displacement {
            result.displacements(2 * node_j) - result.displacements(2 * node_i),
            result.displacements(2 * node_j + 1) -
                result.displacements(2 * node_i + 1)};
        const auto axial_extension =
            dot(state.element_directions[element_index], relative_displacement);
        result.axial_forces[element_index] *= axial_extension;
        result.energies[element_index] *=
            0.5f * axial_extension * axial_extension;
    }

    return result;
}

void make_triangulation(Analysis_state &state)
{
    CDT::Triangulation<float> cdt;

    // FIXME: `vertices` is identical to `nodes`, just a different type
    std::vector<CDT::V2d<float>> vertices;
    vertices.reserve(state.nodes.size());
    for (const auto [x, y] : state.nodes)
    {
        vertices.emplace_back(x, y);
    }
    cdt.insertVertices(vertices);
    cdt.eraseSuperTriangle();

    // FIXME: `edges` is identical to `elements`, just a different type
    const auto edges = CDT::extractEdgesFromTriangles(cdt.triangles);
    state.elements.reserve(edges.size());
    for (const auto &edge : edges)
    {
        const auto [i, j] = edge.verts();
        state.elements.emplace_back(i, j);
    }
}