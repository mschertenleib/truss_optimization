#include "optimization.hpp"

#include "CDT.h"

#include <Eigen/SparseCholesky>

#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>

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

void assemble(Optimization_state &state)
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

void solve_equilibrium_system(Optimization_state &state)
{
    // TODO: test different solvers
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;
    solver.compute(state.stiffness_matrix);
    solver.factorize(state.stiffness_matrix);

    if (const auto result = solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Decomposition failed: {}", to_string(result)));
    }

    const auto num_dofs = state.nodes.size() * 2;
    const auto num_free_dofs = num_dofs - state.fixed_dofs.size();
    Eigen::VectorXf free_displacements(num_free_dofs);
    free_displacements = solver.solve(state.loads);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Solving failed: {}", to_string(result)));
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

void make_triangulation(Optimization_state &state)
{
    CDT::Triangulation<float> cdt;

    cdt.insertVertices(
        state.nodes.cbegin(),
        state.nodes.cend(),
        [](const vec2 &v) { return v.x; },
        [](const vec2 &v) { return v.y; });
    cdt.eraseSuperTriangle();

    // FIXME: is there a better way to do this?
    const auto edges = CDT::extractEdgesFromTriangles(cdt.triangles);
    state.elements.clear();
    state.elements.reserve(edges.size());
    for (const auto &edge : edges)
    {
        const auto [i, j] = edge.verts();
        state.elements.emplace_back(i, j);
    }
}

void equal_stress_projection(Optimization_state &state)
{
    const auto factor =
        max_length /
        state.lengths.cwiseProduct(state.axial_forces.cwiseAbs()).sum();
    state.activations = (state.axial_forces.cwiseAbs() * factor)
                            .cwiseMax(min_activation)
                            .cwiseMin(1.0f);
}

void geometry_step(Optimization_state &state)
{
    // FIXME: this has to be for free DOFs only
    std::vector<vec2> gradients(state.nodes.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto force = state.axial_forces[static_cast<Eigen::Index>(e)];
        const auto gradient_contrib =
            force * force /
            (young_modulus * area *
             state.activations[static_cast<Eigen::Index>(e)]) *
            state.element_directions[e];
        const auto [i, j] = state.elements[e];
        gradients[i] -= gradient_contrib;
        gradients[j] += gradient_contrib;
    }

    auto trial_positions = state.nodes;
    constexpr float gamma {10000.0f};
    constexpr float move_limit {0.02f};
    constexpr unsigned int max_tries {10};

    constexpr auto clamp_to_domain = [](const vec2 &pos)
    {
        // FIXME: this shouldn't be hardcoded
        return vec2 {std::clamp(pos.x, -0.8f, 0.8f),
                     std::clamp(pos.y, -0.8f, 0.8f)};
    };

    // FIXME: don't recompute this
    const auto old_compliance =
        state.loads.dot(state.displacements(state.free_dofs));

    // FIXME: this assumes the first nodes are the fixed ones, followed by
    // the loaded one
    for (std::size_t i {state.fixed_dofs.size() / 2 + 1};
         i < state.nodes.size();
         ++i)
    {
        auto step = -gamma * gradients[i];
        if (const auto norm_step = norm(step); norm_step > move_limit)
        {
            step *= move_limit / norm_step;
        }

        trial_positions[i] = clamp_to_domain(state.nodes[i] + step);
    }

    const auto old_nodes = state.nodes;
    state.nodes = trial_positions;
    assemble(state);
    try
    {
        solve_equilibrium_system(state);
    }
    catch (const std::exception &e)
    {
        state.nodes = old_nodes;
    }

    // TODO: this check might not be necessary?
    const auto new_compliance =
        state.loads.dot(state.displacements(state.free_dofs));
    if (new_compliance >= old_compliance)
    {
        state.nodes = old_nodes;
    }
}

} // namespace

void optimization_init(const std::vector<vec2> &fixed_nodes,
                       const vec2 &load_node,
                       const vec2 &load_vector,
                       Optimization_state &state)
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

    state.fixed_dofs.resize(fixed_nodes.size() * 2);
    std::iota(state.fixed_dofs.begin(), state.fixed_dofs.end(), 0);

    state.free_dofs.resize(2 + num_free_nodes * 2);
    std::iota(state.free_dofs.begin(),
              state.free_dofs.end(),
              state.fixed_dofs.size());

    state.loads.setZero(static_cast<Eigen::Index>(state.free_dofs.size()));
    state.loads(0) = load_vector.x;
    state.loads(1) = load_vector.y;

    state.activations.resize(state.elements.size());
    float total_length {0.0f};
    for (std::size_t element_index {0}; element_index < state.elements.size();
         ++element_index)
    {
        const auto [i, j] = state.elements[element_index];
        total_length += norm(state.nodes[j] - state.nodes[i]);
    }
    std::fill(state.activations.begin(),
              state.activations.end(),
              max_length / total_length);

    assemble(state);
    solve_equilibrium_system(state);
}

void optimization_step(Optimization_state &state)
{
    // Size edges
    equal_stress_projection(state);

    // Move nodes
    geometry_step(state);

    // Add/remove nodes and re-triangulate

    // Solve linear elasticity equilibrium system
    // Currently done in geometry_step, but we might want to only partly compute
    // it there?

    // assemble(state);
    // solve_equilibrium_system(state);
}
