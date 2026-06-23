#include "optimization.hpp"

#include "CDT.h"

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <ranges>

namespace
{

constexpr float young_modulus {200e9f};
constexpr float area {0.000025f};
constexpr float min_activation {1e-3f};
constexpr float max_length {8.0f};

constexpr auto invalid_index = std::numeric_limits<std::uint32_t>::max();

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

void rebuild_sparsity(Optimization_state &state)
{
    state.element_k_indices.resize(state.elements.size());
    state.element_lengths.resize(state.elements.size());
    state.element_directions.resize(state.elements.size());

    std::vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(state.elements.size() * 16);

    std::vector<std::array<std::uint32_t, 4>> element_dofs;
    element_dofs.reserve(state.elements.size());

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto [node_i, node_j] = state.elements[e];
        const std::array<std::uint32_t, 4> e_dofs {
            state.all_to_free_dofs[2 * node_i + 0],
            state.all_to_free_dofs[2 * node_i + 1],
            state.all_to_free_dofs[2 * node_j + 0],
            state.all_to_free_dofs[2 * node_j + 1]};
        element_dofs.push_back(e_dofs);

        for (std::size_t a {0}; a < 4; ++a)
        {
            if (e_dofs[a] == invalid_index)
            {
                continue;
            }

            for (std::size_t b {0}; b < 4; ++b)
            {
                if (e_dofs[b] == invalid_index)
                {
                    continue;
                }

                // NOTE: the value is irrelevant, only the sparsity matters
                triplets.emplace_back(e_dofs[a], e_dofs[b], 1.0f);
            }
        }
    }

    state.stiffness_matrix.resize(
        static_cast<Eigen::Index>(state.free_dofs.size()),
        static_cast<Eigen::Index>(state.free_dofs.size()));
    state.stiffness_matrix.setFromTriplets(triplets.begin(), triplets.end());
    state.stiffness_matrix.makeCompressed();

    // Make the inner indices sorted so we can binary-search them
    // TODO: if our connectivity is very sparse, it might be faster to not sort
    // and do a linear search instead
    state.stiffness_matrix.sortInnerIndices();

    const auto *const outer = state.stiffness_matrix.outerIndexPtr();
    const auto *const inner = state.stiffness_matrix.innerIndexPtr();

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto &e_dofs = element_dofs[e];
        auto &k_indices = state.element_k_indices[e];

        for (std::size_t a {0}; a < 4; ++a)
        {
            for (std::size_t b {0}; b < 4; ++b)
            {
                const auto row = e_dofs[a];
                const auto col = e_dofs[b];
                const auto local_flat_index = 4 * a + b;

                if (row == invalid_index || col == invalid_index)
                {
                    k_indices[local_flat_index] = invalid_index;
                    continue;
                }

                // ColMajor: search within the column's inner-index range
                const auto begin = outer[col];
                const auto end = outer[col + 1];
                const auto *const it =
                    std::lower_bound(inner + begin, inner + end, row);
                assert(it != inner + end &&
                       *it == static_cast<Eigen::Index>(row));

                k_indices[local_flat_index] =
                    static_cast<std::uint32_t>(it - inner);
            }
        }
    }
}

void assemble_stiffness_matrix(Optimization_state &state)
{
    assert(state.elements.size() == state.element_lengths.size());
    assert(state.elements.size() == state.element_directions.size());
    assert(state.elements.size() == state.activations.size());

    auto *const values = state.stiffness_matrix.valuePtr();
    std::fill_n(values, state.stiffness_matrix.nonZeros(), 0.0f);

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto [node_i, node_j] = state.elements[e];
        const auto vec = state.nodes[node_j] - state.nodes[node_i];
        const auto length = norm(vec);
        const auto dir = vec / length;
        state.element_lengths[e] = length;
        state.element_directions[e] = dir;

        const auto stiffness_constant =
            state.activations[e] * young_modulus * area / length;

        const auto [c, s] = dir;
        const auto cc = c * c;
        const auto ss = s * s;
        const auto cs = c * s;
        const float local[16] {cc,
                               cs,
                               -cc,
                               -cs,
                               cs,
                               ss,
                               -cs,
                               -ss,
                               -cc,
                               -cs,
                               cc,
                               cs,
                               -cs,
                               -ss,
                               cs,
                               ss};

        const auto &k_indices = state.element_k_indices[e];
        for (std::size_t l {0}; l < 16; ++l)
        {
            const auto pos = k_indices[l];
            if (pos != invalid_index)
            {
                values[pos] += stiffness_constant * local[l];
            }
        }
    }
}

void solve_equilibrium_system(Optimization_state &state)
{
    state.solver.analyzePattern(state.stiffness_matrix);
    if (const auto result = state.solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(std::format(
            "Symbolic decomposition failed: {}", to_string(result)));
    }

    state.solver.factorize(state.stiffness_matrix);
    if (const auto result = state.solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Numeric decomposition failed: {}", to_string(result)));
    }

    Eigen::VectorXf free_displacements(state.free_dofs.size());
    free_displacements = state.solver.solve(state.loads);
    if (const auto result = state.solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Solving failed: {}", to_string(result)));
    }
    state.displacements.setZero(
        static_cast<Eigen::Index>(state.nodes.size() * 2));
    state.displacements(state.free_dofs) = free_displacements;

    state.axial_forces.resize(state.elements.size());
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
    rebuild_sparsity(state);
    assemble_stiffness_matrix(state);
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

#if 0
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
#else

    constexpr unsigned int num_x {20};
    constexpr unsigned int num_y {20};
    nodes.reserve(num_x * num_y);
    for (unsigned int iy {0}; iy < num_y; ++iy)
    {
        for (unsigned int ix {0}; ix < num_x; ++ix)
        {
            const auto x = -0.8f + static_cast<float>(ix) / (num_x - 1) * 1.6f;
            const auto y = -0.8f + static_cast<float>(iy) / (num_y - 1) * 1.6f;
            nodes.emplace_back(x, y);
        }
    }

    const auto node = [num_x](std::uint32_t iy, std::uint32_t ix)
    { return iy * num_x + ix; };

    const auto dof =
        [&node](std::uint32_t iy, std::uint32_t ix, std::uint32_t axis)
    { return node(iy, ix) * 2 + axis; };

    const std::vector<std::uint32_t> fixed_dofs {dof(num_y / 4, 0, 0),
                                                 dof(num_y / 4, 0, 1),
                                                 dof(num_y * 3 / 4, 0, 0),
                                                 dof(num_y * 3 / 4, 0, 1)};
    const std::vector<std::uint32_t> loaded_dofs {dof(num_y / 2, num_x - 1, 0),
                                                  dof(num_y / 2, num_x - 1, 1)};
    const std::vector<float> loads {0.0f, -1.0f};
    std::vector<std::uint32_t> immovable_dofs;
    immovable_dofs.insert(
        immovable_dofs.end(), fixed_dofs.begin(), fixed_dofs.end());
    immovable_dofs.insert(
        immovable_dofs.end(), loaded_dofs.begin(), loaded_dofs.end());
    std::ranges::sort(immovable_dofs);

    std::vector<Element> elements;
    for (std::uint32_t iy {0}; iy < num_y; ++iy)
    {
        for (std::uint32_t ix {0}; ix < num_x; ++ix)
        {
            if (ix + 1 < num_x)
                elements.emplace_back(node(iy, ix), node(iy, ix + 1));
            if (iy + 1 < num_y)
                elements.emplace_back(node(iy, ix), node(iy + 1, ix));
            if ((ix + 1 < num_x) && (iy + 1 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 1, ix + 1));
            if ((ix > 0) && (iy + 1 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 1, ix - 1));
            if ((ix + 2 < num_x) && (iy + 1 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 1, ix + 2));
            if ((ix + 1 < num_x) && (iy + 2 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 2, ix + 1));
            if ((ix > 0) && (iy + 2 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 2, ix - 1));
            if ((ix > 1) && (iy + 1 < num_y))
                elements.emplace_back(node(iy, ix), node(iy + 1, ix - 2));
        }
    }

#endif

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
    // FIXME: this must change place, and we need to rebuild sparsity if
    // topology changes
    assemble_stiffness_matrix(state);
    solve_equilibrium_system(state);
}
