#include "optimization.hpp"

#include "CDT.h"

#include <Eigen/SparseCholesky>

#include <algorithm>
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

constexpr auto invalid_index = static_cast<std::uint32_t>(-1);

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
    assert(state.elements.size() == state.activations.size());

    std::vector<Eigen::Triplet<float>> triplets;

    state.element_lengths.resize(state.elements.size());
    state.element_directions.resize(state.elements.size());

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto [node_i, node_j] = state.elements[e];
        const auto activation = state.activations[e];
        const auto vec = state.nodes[node_j] - state.nodes[node_i];
        const auto length = norm(vec);
        state.element_lengths[e] = length;
        const auto dir = vec / length;
        state.element_directions[e] = dir;
        const auto stiffness_constant =
            activation * (young_modulus * area) / length;

        const auto [c, s] = dir;
        const float element_stiffness_matrix[4][4] {
            {c * c, c * s, -c * c, -c * s},
            {c * s, s * s, -c * s, -s * s},
            {-c * c, -c * s, c * c, c * s},
            {-c * s, -s * s, c * s, s * s}};

        const std::uint32_t dof_indices[4] {
            state.all_to_free_dofs[2 * node_i],
            state.all_to_free_dofs[2 * node_i + 1],
            state.all_to_free_dofs[2 * node_j],
            state.all_to_free_dofs[2 * node_j + 1]};
        for (unsigned int a {0}; a < 4; ++a)
        {
            for (unsigned int b {0}; b < 4; ++b)
            {
                if (dof_indices[a] != invalid_index &&
                    dof_indices[b] != invalid_index)
                {
                    triplets.emplace_back(dof_indices[a],
                                          dof_indices[b],
                                          stiffness_constant *
                                              element_stiffness_matrix[a][b]);
                }
            }
        }
    }

    state.stiffness_matrix.setZero();
    state.stiffness_matrix.resize(
        static_cast<Eigen::Index>(state.free_dofs.size()),
        static_cast<Eigen::Index>(state.free_dofs.size()));
    state.stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    state.stiffness_matrix.prune(0.0f);
}

void solve_equilibrium_system(Optimization_state &state)
{
    // TODO: test different solvers
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;

    solver.analyzePattern(state.stiffness_matrix);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(std::format(
            "Symbolic decomposition failed: {}", to_string(result)));
    }

    solver.factorize(state.stiffness_matrix);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Numeric decomposition failed: {}", to_string(result)));
    }

    Eigen::VectorXf free_displacements(state.free_dofs.size());
    free_displacements = solver.solve(state.loads);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Solving failed: {}", to_string(result)));
    }
    state.displacements.setZero(
        static_cast<Eigen::Index>(state.nodes.size() * 2));
    state.displacements(state.free_dofs) = free_displacements;

    state.axial_forces.resize(state.elements.size());
    state.energies.resize(state.elements.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto [node_i, node_j] = state.elements[e];
        const vec2 relative_displacement {
            state.displacements(2 * node_j) - state.displacements(2 * node_i),
            state.displacements(2 * node_j + 1) -
                state.displacements(2 * node_i + 1)};
        const auto axial_extension =
            dot(state.element_directions[e], relative_displacement);
        const auto stiffness_constant = state.activations[e] *
                                        (young_modulus * area) /
                                        state.element_lengths[e];

        state.axial_forces[e] = stiffness_constant * axial_extension;
        state.energies[e] =
            0.5f * stiffness_constant * axial_extension * axial_extension;
    }
}

void make_triangulation(const std::vector<vec2> &nodes,
                        std::vector<Element> &elements)
{
    CDT::Triangulation<float> cdt;
    cdt.insertVertices(
        nodes.cbegin(),
        nodes.cend(),
        [](const vec2 &v) { return v.x; },
        [](const vec2 &v) { return v.y; });
    cdt.eraseSuperTriangle();

    // FIXME: is there a better way to do this?
    const auto edges = CDT::extractEdgesFromTriangles(cdt.triangles);
    elements.clear();
    elements.reserve(edges.size());
    for (const auto &edge : edges)
    {
        const auto [i, j] = edge.verts();
        elements.emplace_back(i, j);
    }
}

void optimization_init(const std::vector<vec2> &nodes,
                       const std::vector<Element> &elements,
                       const std::vector<std::uint32_t> &fixed_dofs,
                       const std::vector<std::uint32_t> &loaded_dofs,
                       const std::vector<float> &loads,
                       const std::vector<std::uint32_t> &immovable_dofs,
                       Optimization_state &state)
{
    const auto num_dofs = nodes.size() * 2;

    assert(std::ranges::all_of(elements,
                               [num_nodes = nodes.size()](auto el)
                               {
                                   return el.node_i < num_nodes &&
                                          el.node_j < num_nodes &&
                                          el.node_i != el.node_j;
                               }));

    assert(std::ranges::is_sorted(fixed_dofs));
    assert(std::ranges::adjacent_find(fixed_dofs) ==
           std::ranges::end(fixed_dofs));
    assert(fixed_dofs.empty() || fixed_dofs.back() < num_dofs);

    assert(loaded_dofs.size() == loads.size());
    assert(std::ranges::all_of(
        loaded_dofs, [num_dofs](auto dof) { return dof < num_dofs; }));

    assert(std::ranges::is_sorted(immovable_dofs));
    assert(std::ranges::adjacent_find(immovable_dofs) ==
           std::ranges::end(immovable_dofs));
    assert(immovable_dofs.empty() || immovable_dofs.back() < num_dofs);

    // FIXME: avoid re-allocations

    state.nodes = nodes;
    state.elements = elements;

    state.fixed_dofs = fixed_dofs;
    state.immovable_dofs = immovable_dofs;

    const auto num_free_dofs = num_dofs - fixed_dofs.size();
    state.free_dofs.resize(num_free_dofs);
    state.all_to_free_dofs.resize(num_dofs);
    std::uint32_t free_dof_index {0};
    std::uint32_t fixed_dof_index {0};
    for (std::uint32_t i {0}; i < num_dofs; ++i)
    {
        if (fixed_dof_index < fixed_dofs.size() &&
            i == fixed_dofs[fixed_dof_index])
        {
            state.all_to_free_dofs[i] = invalid_index;
            ++fixed_dof_index;
        }
        else
        {
            state.all_to_free_dofs[i] = free_dof_index;
            state.free_dofs[free_dof_index] = i;
            ++free_dof_index;
        }
    }

    state.loads.setZero(static_cast<Eigen::Index>(state.free_dofs.size()));
    for (const auto [loaded_dof, load] : std::views::zip(loaded_dofs, loads))
    {
        const auto index = state.all_to_free_dofs[loaded_dof];
        if (index != invalid_index)
        {
            state.loads(index) = load;
        }
    }

    float total_length {0.0f};
    for (const auto [i, j] : state.elements)
    {
        total_length += norm(state.nodes[j] - state.nodes[i]);
    }
    state.activations.resize(state.elements.size());
    std::ranges::fill(state.activations, max_length / total_length);

    // FIXME: this must disappear from here
    assemble(state);
    solve_equilibrium_system(state);
}

void compute_activations(Optimization_state &state)
{
    // Compute activations by equal stress projection, i.e. the activations that
    // would cause the structure to be fully equally stressed.

    assert(state.element_lengths.size() == state.axial_forces.size());
    assert(state.element_lengths.size() == state.activations.size());

    float sum {0.0f};
    for (std::size_t e {0}; e < state.activations.size(); ++e)
    {
        sum += state.element_lengths[e] * std::abs(state.axial_forces[e]);
    }
    const auto factor = max_length / sum;

    for (std::size_t e {0}; e < state.activations.size(); ++e)
    {
        state.activations[e] = std::clamp(
            std::abs(state.axial_forces[e]) * factor, min_activation, 1.0f);
    }
}

void geometry_step(Optimization_state &state)
{
    state.gradients.clear();
    state.gradients.resize(state.nodes.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto t = state.element_directions[e];
        const auto L = state.element_lengths[e];
        const auto [node_i, node_j] = state.elements[e];
        const vec2 u_i {state.displacements(2 * node_i),
                        state.displacements(2 * node_i + 1)};
        const vec2 u_j {state.displacements(2 * node_j),
                        state.displacements(2 * node_j + 1)};
        const auto u_rel = u_j - u_i;
        const auto s = dot(t, u_rel);
        const vec2 P_u_rel {(1.0f - t.x * t.x) * u_rel.x - t.x * t.y * u_rel.y,
                            -t.x * t.y * u_rel.x +
                                (1.0f - t.y * t.y) * u_rel.y};

        const auto gradient_contrib = young_modulus * area *
                                      state.activations[e] / (L * L) *
                                      (t * (s * s) - (2.0f * s) * P_u_rel);

        state.gradients[node_i] -= gradient_contrib;
        state.gradients[node_j] += gradient_contrib;
    }

    constexpr float gamma {10000.0f};
    constexpr float move_limit {0.02f};

    constexpr auto clamp_to_domain = [](const vec2 &pos)
    {
        // FIXME: this shouldn't be hardcoded
        return vec2 {std::clamp(pos.x, -0.8f, 0.8f),
                     std::clamp(pos.y, -0.8f, 0.8f)};
    };

    std::uint32_t immovable_dof_index {0};
    const auto is_next_immovable_dof =
        [&immovable_dof_index, &state](std::uint32_t dof)
    {
        if (immovable_dof_index < state.immovable_dofs.size() &&
            dof == state.immovable_dofs[immovable_dof_index])
        {
            ++immovable_dof_index;
            return true;
        }
        return false;
    };
    for (std::uint32_t i {0}; i < state.nodes.size(); ++i)
    {
        auto step = -gamma * state.gradients[i];
        if (const auto norm_step = norm(step); norm_step > move_limit)
        {
            step *= move_limit / norm_step;
        }

        const auto new_pos = clamp_to_domain(state.nodes[i] + step);
        if (!is_next_immovable_dof(i * 2))
        {
            state.nodes[i].x = new_pos.x;
        }
        if (!is_next_immovable_dof(i * 2 + 1))
        {
            state.nodes[i].y = new_pos.y;
        }
    }
}

} // namespace

void optimization_create_problem(Optimization_state &state)
{
    std::vector<vec2> nodes;

    nodes.emplace_back(-0.8f, -0.4f);
    nodes.emplace_back(-0.8f, 0.4f);
    const std::vector<std::uint32_t> fixed_dofs {0, 1, 2, 3};

    nodes.emplace_back(0.8f, 0.0f);
    const std::vector<std::uint32_t> loaded_dofs {4, 5};
    const std::vector<float> loads {0.0f, -1.0f};

    const std::vector<std::uint32_t> immovable_dofs {0, 1, 2, 3, 4, 5};

    constexpr std::uint32_t num_free_nodes {200};
    std::minstd_rand rng(17657575);
    std::uniform_real_distribution<float> x(-0.8f, 0.8f);
    std::uniform_real_distribution<float> y(-0.8f, 0.8f);
    nodes.reserve(nodes.size() + num_free_nodes);
    for (std::uint32_t i {0}; i < num_free_nodes; ++i)
    {
        nodes.emplace_back(x(rng), y(rng));
    }

    std::vector<Element> elements;
    make_triangulation(nodes, elements);
    /*for (std::uint32_t i {0}; i < nodes.size(); ++i)
    {
        for (std::uint32_t j {i + 1}; j < nodes.size(); ++j)
        {
            elements.emplace_back(i, j);
        }
    }*/

    optimization_init(
        nodes, elements, fixed_dofs, loaded_dofs, loads, immovable_dofs, state);
}

void optimization_step(Optimization_state &state)
{
    /*
    Full optimization procedure:
    - Initial node distribution
    - Triangulate (there might be a way to generate blue noise and a
      triangulation as a single pass, might be worth looking into it)
    - Loop:
        - If topology changed:
            - Assemble K (might be possible to only update coefficients with a
              topology change)
            - Symbolically decompose K
        - Else:
            - Assemble K (only update values, same sparsity)
        - Numerically decompose K
        - Solve linear system
        - Compute activations
        - Move nodes
        - If necessary, re-triangulate
    */

    // Size edges
    compute_activations(state);

    // Move nodes
    geometry_step(state);

    // Add/remove nodes and re-triangulate

    // Solve linear elasticity equilibrium system
    assemble(state);
    solve_equilibrium_system(state);
}
