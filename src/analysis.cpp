#include "analysis.hpp"

#include "CDT.h"

#include <Eigen/SparseCholesky>

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

namespace
{

constexpr float young_modulus {200e9f};
constexpr float area {0.000025f};
constexpr float min_activation {1e-3f};
constexpr float max_length {8.0f};

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

void assemble(Analysis_state &state)
{
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

    state.lengths.resize(state.elements.size());
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
        state.lengths[element_index] = length;
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

void solve_equilibrium_system(Analysis_state &state)
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

    state.displacements.setZero(num_dofs);
    state.displacements(state.free_dofs) = free_displacements;

    state.axial_forces.resize(
        static_cast<Eigen::Index>(state.stiffness_constants.size()));
    state.energies.resize(
        static_cast<Eigen::Index>(state.stiffness_constants.size()));
    for (std::size_t i {0}; i < state.stiffness_constants.size(); ++i)
    {
        state.axial_forces[static_cast<Eigen::Index>(i)] =
            state.stiffness_constants[i];
        state.energies[static_cast<Eigen::Index>(i)] =
            state.stiffness_constants[i];
    }
    for (std::size_t element_index {0}; element_index < state.elements.size();
         ++element_index)
    {
        const auto [node_i, node_j] = state.elements[element_index];
        const vec2 relative_displacement {
            state.displacements(2 * node_j) - state.displacements(2 * node_i),
            state.displacements(2 * node_j + 1) -
                state.displacements(2 * node_i + 1)};
        const auto axial_extension =
            dot(state.element_directions[element_index], relative_displacement);
        state.axial_forces[static_cast<Eigen::Index>(element_index)] *=
            axial_extension;
        state.energies[static_cast<Eigen::Index>(element_index)] *=
            0.5f * axial_extension * axial_extension;
    }
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
    state.elements.clear();
    state.elements.reserve(edges.size());
    for (const auto &edge : edges)
    {
        const auto [i, j] = edge.verts();
        state.elements.emplace_back(i, j);
    }
}

void equal_stress_projection(Analysis_state &state)
{
    const auto equal_stress =
        state.lengths.cwiseProduct(state.axial_forces.cwiseAbs()).sum() *
        (1.0f / (area * max_length));
    state.activations = (state.axial_forces.cwiseAbs() / (equal_stress * area))
                            .cwiseMax(min_activation)
                            .cwiseMin(1.0f);
}

void geometry_step(Analysis_state &state)
{
    // FIXME: this has to be for free nodes only
    std::vector<vec2> gradients(state.nodes.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto gradient_contrib =
            2.0f * state.energies[static_cast<Eigen::Index>(e)] /
            state.lengths[static_cast<Eigen::Index>(e)] *
            state.element_directions[e];
        const auto [i, j] = state.elements[e];
        gradients[i] -= gradient_contrib;
        gradients[j] += gradient_contrib;
    }

    for (const auto [x, y] : gradients)
    {
        std::cout << x << ' ' << y << '\n';
    }
}

} // namespace

void optimization_init(const std::vector<vec2> &fixed_nodes,
                       const vec2 &load_node,
                       const vec2 &load_vector,
                       Analysis_state &state)
{
    state.nodes.clear();
    state.nodes.insert(
        state.nodes.cbegin(), fixed_nodes.cbegin(), fixed_nodes.cend());
    state.nodes.push_back(load_node);

    constexpr std::uint32_t num_free_nodes {100};
    std::minstd_rand rng(17657575);
    std::uniform_real_distribution<float> x(-0.8f, 0.8f);
    std::uniform_real_distribution<float> y(-0.8f, 0.8f);
    for (std::uint32_t i {0}; i < num_free_nodes; ++i)
    {
        state.nodes.push_back({x(rng), y(rng)});
    }

    make_triangulation(state);

    state.activations.resize(state.elements.size());
    std::fill(state.activations.begin(), state.activations.end(), 1.0f);

    state.fixed_dofs.resize(fixed_nodes.size() * 2);
    std::iota(state.fixed_dofs.begin(), state.fixed_dofs.end(), 0);

    state.free_dofs.resize(2 + num_free_nodes * 2);
    std::iota(state.free_dofs.begin(),
              state.free_dofs.end(),
              state.fixed_dofs.size());

    state.loads.setZero(static_cast<Eigen::Index>(state.free_dofs.size()));
    state.loads(0) = load_vector.x;
    state.loads(1) = load_vector.y;

    // Solve linear elasticity equilibrium system
    assemble(state);
    solve_equilibrium_system(state);
}

void optimization_step(Analysis_state &state)
{
    // Sizing step (size edges)
    equal_stress_projection(state);

    // Geometry step (move nodes)
    geometry_step(state);

    // Add/remove nodes and re-triangulate

    // Solve linear elasticity equilibrium system
    assemble(state);
    solve_equilibrium_system(state);
}
