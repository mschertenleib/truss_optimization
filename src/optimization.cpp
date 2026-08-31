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
constexpr float min_stiffness_activation {1e-5f};
constexpr float activation_move_limit {0.25f};
constexpr float penalization {2.0f};
constexpr float max_area {0.02f * 0.02f};
constexpr float max_volume {10.0f * max_area};
constexpr float min_coordinate {-0.8f};
constexpr float max_coordinate {0.8f};

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

void rebuild_stiffness_matrix_sparsity(Optimization_state &state)
{
    const auto num_elements = state.elements.size();

    std::vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(num_elements * 16);

    std::vector<std::array<std::uint32_t, 4>> element_dofs;
    element_dofs.reserve(num_elements);

    for (std::size_t e {0}; e < num_elements; ++e)
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

    state.element_k_indices.resize(num_elements);

    const auto *const outer = state.stiffness_matrix.outerIndexPtr();
    const auto *const inner = state.stiffness_matrix.innerIndexPtr();

    for (std::size_t e {0}; e < num_elements; ++e)
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

    state.element_directions.resize(num_elements);
    state.element_lengths.resize(num_elements);
    state.activations.resize(num_elements);
}

void assemble_stiffness_matrix(Optimization_state &state)
{
    assert(state.elements.size() == state.element_k_indices.size());
    assert(state.elements.size() == state.element_directions.size());
    assert(state.elements.size() == state.element_lengths.size());
    assert(state.elements.size() == state.activations.size());

    auto *const values = state.stiffness_matrix.valuePtr();
    std::fill_n(values, state.stiffness_matrix.nonZeros(), 0.0f);

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto [i, j] = state.elements[e];
        const auto vec = state.nodes[j] - state.nodes[i];
        const auto length = norm(vec);
        const auto dir = vec / length;
        state.element_lengths[e] = length;
        state.element_directions[e] = dir;

        const auto stiffness_interpolation =
            min_stiffness_activation +
            std::pow(state.activations[e], penalization) *
                (1.0f - min_stiffness_activation);
        const auto stiffness_constant =
            stiffness_interpolation * max_area * young_modulus / length;

        const auto [c, s] = dir;
        const auto cc = c * c;
        const auto ss = s * s;
        const auto cs = c * s;
        const std::array<float, 16> local {cc,
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
        for (std::size_t k {0}; k < 16; ++k)
        {
            const auto pos = k_indices[k];
            if (pos != invalid_index)
            {
                values[pos] += stiffness_constant * local[k];
            }
        }
    }
}

void symbolically_decompose_stiffness_matrix(Optimization_state &state)
{
    state.solver.analyzePattern(state.stiffness_matrix);
    if (const auto result = state.solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(std::format(
            "Symbolic decomposition failed: {}", to_string(result)));
    }
}

void numerically_decompose_stiffness_matrix(Optimization_state &state)
{
    state.solver.factorize(state.stiffness_matrix);
    if (const auto result = state.solver.info(); result != Eigen::Success)
    {
        throw std::runtime_error(
            std::format("Numeric decomposition failed: {}", to_string(result)));
    }
}

void solve_equilibrium_system(Optimization_state &state)
{
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
}

// FIXME: change function name
void compute_compliance_and_gradients(Optimization_state &state)
{
    // FIXME: move all displacement-derived values here (axial_forces, deltas if
    // necessary, etc.)

    state.axial_forces.resize(
        state.elements.size()); // TODO: these should probably be sized only
                                // when elements (connectivity) change
    state.compliance = 0.0f;
    state.dC_dx.clear(); // TODO: this should really just be a fill with zeros
    state.dC_dx.resize(state.nodes.size());

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto dir = state.element_directions[e];
        const auto length = state.element_lengths[e];
        const auto [i, j] = state.elements[e];
        const vec2 u_i {state.displacements(2 * i),
                        state.displacements(2 * i + 1)};
        const vec2 u_j {state.displacements(2 * j),
                        state.displacements(2 * j + 1)};
        const auto area = state.activations[e] * max_area;
        const auto delta_u = u_j - u_i;
        const auto delta = dot(dir, delta_u);

        const auto force = area * young_modulus / length * delta;
        state.axial_forces[e] = force;
        state.compliance += force * delta;

        const auto delta_u_transverse = delta_u - delta * dir;
        const auto dCe_dxj =
            area * young_modulus / (length * length) *
            ((delta * delta) * dir - (2.0f * delta) * delta_u_transverse);
        state.dC_dx[j] += dCe_dxj;
        state.dC_dx[i] -= dCe_dxj;
    }
}

void compute_areas(Optimization_state &state)
{
    const auto sqrt_factor = std::sqrt((1.0f - min_stiffness_activation) *
                                       penalization * young_modulus);
    std::vector<float> r(state.elements.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto dir = state.element_directions[e];
        const auto length = state.element_lengths[e];
        const auto [i, j] = state.elements[e];
        const vec2 u_i {state.displacements(2 * i),
                        state.displacements(2 * i + 1)};
        const vec2 u_j {state.displacements(2 * j),
                        state.displacements(2 * j + 1)};
        const auto delta_u = u_j - u_i;
        const auto delta = dot(dir, delta_u);

        // FIXME: this might be wrong with the new stiffness interpolation.
        // Currently it is 0 if activations are 0, but there must be a way to
        // "re-activate" activations that were 0. Also, we should check if this
        // is also a problem in the MATLAB code, and how they do it, since it
        // ends up being very similar to what we do.
        r[e] = sqrt_factor * std::abs(delta) / length *
               std::pow(state.activations[e], (penalization + 1.0f) * 0.5f);
    }

    const auto compute_volume = [&](float mu)
    {
        float volume {0.0f};

        for (std::size_t e {0}; e < state.elements.size(); ++e)
        {
            const auto lower =
                state.activations[e] / (1.0f + activation_move_limit);
            const auto upper =
                state.activations[e] * (1.0f + activation_move_limit);
            const auto new_activation =
                std::clamp(std::clamp(r[e] / mu, lower, upper), 0.0f, 1.0f);
            volume += new_activation * max_area * state.element_lengths[e];
        }

        return volume;
    };

    float mu_lo {0.0f};
    float mu_hi {1.0f};

    // Increase mu until the volume is below the target
    int bound_iter {0};
    while (compute_volume(mu_hi) > max_volume)
    {
        mu_hi *= 2.0f;
        ++bound_iter;
        if (bound_iter >= 20)
        {
            std::cerr << "Volume constraint not feasible\n";
            break;
        }
    }

    for (int iter {0}; iter < 100; ++iter)
    {
        const auto mu_mid = 0.5f * (mu_lo + mu_hi);

        if (compute_volume(mu_mid) > max_volume)
        {
            mu_lo = mu_mid;
        }
        else
        {
            mu_hi = mu_mid;
        }

        if (mu_hi - mu_lo <= 1e-4f * (1.0f + mu_hi))
        {
            break;
        }
    }
    const auto mu = 0.5f * (mu_lo + mu_hi);
    std::cout << "mu=" << mu << '\n';

    state.new_activations.resize(state.elements.size());
    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        // FIXME: code duplication
        const auto lower =
            state.activations[e] / (1.0f + activation_move_limit);
        const auto upper =
            state.activations[e] * (1.0f + activation_move_limit);
        state.new_activations[e] =
            std::clamp(std::clamp(r[e] / mu, lower, upper), 0.0f, 1.0f);
    }
}

void update_geometry(Optimization_state &state)
{
    float q_max {0.0f};
    for (const auto &g : state.dC_dx)
    {
        const auto q = -g;
        q_max = std::max(norm(q), q_max);
    }

    state.predicted_compliance_reduction = 0.0f;
    for (const auto i : state.movable_nodes)
    {
        const auto x = state.nodes[i];
        const auto g = state.dC_dx[i];
        const auto q = -g;
        const auto q_hat = q / q_max;
        const auto x_raw = x + state.delta_x * q_hat;

        vec2 x_cand {std::clamp(x_raw.x, min_coordinate, max_coordinate),
                     std::clamp(x_raw.y, min_coordinate, max_coordinate)};

        state.predicted_compliance_reduction -= dot(g, x_cand - x);
        state.nodes[i] = x_cand;
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
                       const std::vector<std::uint32_t> &immovable_nodes,
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

    assert(std::ranges::is_sorted(immovable_nodes));
    assert(std::ranges::adjacent_find(immovable_nodes) ==
           std::ranges::end(immovable_nodes));
    assert(immovable_nodes.empty() || immovable_nodes.back() < nodes.size());

    // FIXME: avoid re-allocations

    state.nodes = nodes;
    state.elements = elements;

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

    state.movable_nodes.resize(nodes.size() - immovable_nodes.size());
    std::uint32_t movable_node_index {0};
    std::uint32_t immovable_node_index {0};
    for (std::uint32_t i {0}; i < nodes.size(); ++i)
    {
        if (immovable_node_index < immovable_nodes.size() &&
            i == immovable_nodes[immovable_node_index])
        {
            ++immovable_node_index;
        }
        else
        {
            state.movable_nodes[movable_node_index] = i;
            ++movable_node_index;
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
    std::ranges::fill(state.activations,
                      max_volume / (max_area * total_length));

    state.sparsity_stale = true;

    const auto mean_length =
        total_length / static_cast<float>(state.elements.size());
    state.delta_x = 0.05f * mean_length;
    state.step = 0;
}

} // namespace

void optimization_create_problem(Optimization_state &state, Problem problem)
{
    if (problem == Problem::random_delaunay)
    {
        std::vector<vec2> nodes;
        nodes.emplace_back(-0.8f, -0.4f);
        nodes.emplace_back(-0.8f, 0.4f);
        const std::vector<std::uint32_t> fixed_dofs {0, 1, 2, 3};

        nodes.emplace_back(0.8f, 0.0f);
        const std::vector<std::uint32_t> loaded_dofs {4, 5};
        const std::vector<float> loads {0.0f, -1.0f};

        const std::vector<std::uint32_t> immovable_nodes {0, 1, 2};

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

        optimization_init(nodes,
                          elements,
                          fixed_dofs,
                          loaded_dofs,
                          loads,
                          immovable_nodes,
                          state);
    }
    else if (problem == Problem::regular_grid)
    {
        std::vector<vec2> nodes;
        constexpr unsigned int num_x {20};
        constexpr unsigned int num_y {20};
        nodes.reserve(num_x * num_y);
        for (unsigned int iy {0}; iy < num_y; ++iy)
        {
            for (unsigned int ix {0}; ix < num_x; ++ix)
            {
                const auto x =
                    -0.8f + static_cast<float>(ix) / (num_x - 1) * 1.6f;
                const auto y =
                    -0.8f + static_cast<float>(iy) / (num_y - 1) * 1.6f;
                nodes.emplace_back(x, y);
            }
        }

        const auto node = [num_x](std::uint32_t iy, std::uint32_t ix)
        { return iy * num_x + ix; };

        const std::array fixed_nodes {node(num_y / 4, 0),
                                      node(num_y * 3 / 4, 0)};
        const std::vector<std::uint32_t> fixed_dofs {fixed_nodes[0] * 2 + 0,
                                                     fixed_nodes[0] * 2 + 1,
                                                     fixed_nodes[1] * 2 + 0,
                                                     fixed_nodes[1] * 2 + 1};

        const auto loaded_node = node(num_y / 2, num_x - 1);
        const std::vector<std::uint32_t> loaded_dofs {loaded_node * 2 + 1};
        const std::vector<float> loads {-1.0f};

        std::vector<std::uint32_t> immovable_nodes {
            fixed_nodes[0], fixed_nodes[1], loaded_node};
        std::ranges::sort(immovable_nodes);

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

        optimization_init(nodes,
                          elements,
                          fixed_dofs,
                          loaded_dofs,
                          loads,
                          immovable_nodes,
                          state);
    }
    else if (problem == Problem::hexagonal)
    {
        constexpr unsigned int num_x {20};
        const auto side = 1.6f / static_cast<float>(num_x - 1);
        const auto height = std::numbers::sqrt3_v<float> * 0.5f * side;
        const auto num_y =
            static_cast<unsigned int>(std::floor(1.6f / height)) + 1;
        const auto total_height = static_cast<float>(num_y - 1) * height;

        std::vector<vec2> nodes;
        for (unsigned int iy {0}; iy < num_y; ++iy)
        {
            const auto y =
                -total_height * 0.5f + static_cast<float>(iy) * height;
            const auto is_offset = static_cast<bool>(iy & 1);
            for (unsigned int ix {0}; ix < (is_offset ? num_x - 1 : num_x);
                 ++ix)
            {
                const auto x = -0.8f + (is_offset ? 0.5f * side : 0.0f) +
                               static_cast<float>(ix) * side;

                nodes.emplace_back(x, y);
            }
        }

        const auto closest = [&nodes](float x, float y)
        {
            float min_distance {std::numeric_limits<float>::max()};
            std::uint32_t index {};
            for (std::uint32_t i {0}; i < nodes.size(); ++i)
            {
                if (const auto distance = norm(vec2 {x, y} - nodes[i]);
                    distance < min_distance)
                {
                    min_distance = distance;
                    index = i;
                }
            }
            return index;
        };

        const std::array fixed_nodes {closest(-0.8f, -0.4f),
                                      closest(-0.8f, 0.4f)};
        const std::vector<std::uint32_t> fixed_dofs {fixed_nodes[0] * 2 + 0,
                                                     fixed_nodes[0] * 2 + 1,
                                                     fixed_nodes[1] * 2 + 0,
                                                     fixed_nodes[1] * 2 + 1};

        const auto loaded_node = closest(0.8f, 0.0f);
        const std::vector<std::uint32_t> loaded_dofs {loaded_node * 2 + 1};
        const std::vector<float> loads {-1.0f};

        std::vector<std::uint32_t> immovable_nodes {
            fixed_nodes[0], fixed_nodes[1], loaded_node};
        std::ranges::sort(immovable_nodes);

        std::vector<Element> elements;
        make_triangulation(nodes, elements);

        optimization_init(nodes,
                          elements,
                          fixed_dofs,
                          loaded_dofs,
                          loads,
                          immovable_nodes,
                          state);
    }
    else
    {
        std::cerr << std::format("Unhandled problem: value = {}\n",
                                 std::to_underlying(problem));
    }
}

void optimization_step(Optimization_state &state)
{
    if (state.sparsity_stale)
    {
        rebuild_stiffness_matrix_sparsity(state);
        symbolically_decompose_stiffness_matrix(state);
        state.sparsity_stale = false;
    }

    assemble_stiffness_matrix(state);
    numerically_decompose_stiffness_matrix(state);
    solve_equilibrium_system(state);

    const auto prev_compliance = state.compliance;
    compute_compliance_and_gradients(state);
    if (state.step > 0)
    {
        const auto compliance_reduction = prev_compliance - state.compliance;
        const auto trust_region_ratio =
            compliance_reduction / state.predicted_compliance_reduction;
        constexpr float r_low {0.25f};
        constexpr float r_high {0.75f};
        constexpr float gamma_dec {0.5f};
        constexpr float gamma_inc {1.5f};
        if (trust_region_ratio < r_low)
            state.delta_x = gamma_dec * state.delta_x;
        else if (trust_region_ratio > r_high)
            state.delta_x = gamma_inc * state.delta_x;

        std::cout << std::format(
            "ratio={} delta={}\n", trust_region_ratio, state.delta_x);
    }

    compute_areas(state);

    // update_geometry(state);

    float max_diff {0.0f};
    for (std::uint32_t i {0}; i < state.nodes.size(); ++i)
    {
        const auto diff =
            std::abs(state.new_activations[i] - state.activations[i]);
        if (diff > max_diff)
        {
            max_diff = diff;
        }
    }
    std::cout << "|A_new-A|=" << max_diff << '\n';

    state.activations = state.new_activations;

    ++state.step;
}
