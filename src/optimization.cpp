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
constexpr float max_area {0.005f * 0.005f};
constexpr float min_area {1e-3f * max_area};
constexpr float max_volume {8.0f * max_area};
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
    assert(state.elements.size() == state.areas.size());

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

        const auto stiffness_constant = state.areas[e] * young_modulus / length;

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

void compute_gradients(Optimization_state &state)
{
    state.compliance = 0.0f;
    state.volume = 0.0f;

    state.dC_dA.resize(state.elements.size());
    state.dC_dx.clear();
    state.dC_dx.resize(state.nodes.size());

    state.dV_dA.resize(state.elements.size());
    state.dV_dx.clear();
    state.dV_dx.resize(state.nodes.size());

    for (std::size_t e {0}; e < state.elements.size(); ++e)
    {
        const auto dir = state.element_directions[e];
        const auto length = state.element_lengths[e];
        const auto [i, j] = state.elements[e];
        const vec2 u_i {state.displacements(2 * i),
                        state.displacements(2 * i + 1)};
        const vec2 u_j {state.displacements(2 * j),
                        state.displacements(2 * j + 1)};
        const auto area = state.areas[e];
        const auto delta_u = u_j - u_i;
        const auto delta = dot(dir, delta_u);

        state.compliance += area * young_modulus / length * delta * delta;
        state.volume += area * length;

        state.dC_dA[e] = -young_modulus / length * delta * delta;

        const auto delta_u_transverse = delta_u - delta * dir;
        const auto dCe_dxj =
            area * young_modulus / (length * length) *
            ((delta * delta) * dir - (2.0f * delta) * delta_u_transverse);

        state.dC_dx[j] += dCe_dxj;
        state.dC_dx[i] -= dCe_dxj;

        state.dV_dA[e] = length;

        const auto dVe_dxj = area * dir;
        state.dV_dx[j] += dVe_dxj;
        state.dV_dx[i] -= dVe_dxj;
    }
}

void mma_init(MMA_state &mma,
              const Eigen::ArrayXf &x_min,
              const Eigen::ArrayXf &x_max,
              const MMA_options &options)
{
    assert(x_min.size() == x_max.size());
    assert((x_min <= x_max).all());

    mma.size = x_min.size();

    mma.x_min = x_min;
    mma.x_max = x_max;
    mma.range = x_max - x_min;
    mma.options = options;

    mma.s_old_1.resize(mma.size);
    mma.s_old_2.resize(mma.size);
    mma.lower.resize(mma.size);
    mma.upper.resize(mma.size);

    mma.iteration = 0;
    mma.initialized = false;
}

void mma_update_asymptotes(MMA_state &mma, const Eigen::ArrayXf &s)
{
    for (Eigen::Index i {0}; i < mma.size; ++i)
    {
        const auto product =
            (s[i] - mma.s_old_1[i]) * (mma.s_old_1[i] - mma.s_old_2[i]);

        float factor {1.0f};
        if (product > 0.0f)
        {
            factor = mma.options.asymptote_inc;
        }
        else if (product < 0.0f)
        {
            factor = mma.options.asymptote_dec;
        }

        mma.lower[i] = s[i] - factor * (mma.s_old_1[i] - mma.lower[i]);
        mma.upper[i] = s[i] + factor * (mma.upper[i] - mma.s_old_1[i]);

        const auto lower_min = s[i] - 10.0f;
        const auto lower_max = s[i] - 0.01f;
        const auto upper_min = s[i] + 0.01f;
        const auto upper_max = s[i] + 10.0f;
        mma.lower[i] = std::max(mma.lower[i], lower_min);
        mma.lower[i] = std::min(mma.lower[i], lower_max);
        mma.upper[i] = std::max(mma.upper[i], upper_min);
        mma.upper[i] = std::min(mma.upper[i], upper_max);
    }
}

/*
 * One MMA outer iteration.
 *
 * Input:
 *   x      : current physical design
 *   f      : objective value f(x)
 *   df     : physical gradient df/dx
 *   g      : scalar constraint, g(x) <= 0
 *   dg     : physical gradient dg/dx
 *
 * Output:
 *   x_new  : next physical design
 *
 * Requirement:
 *   current design should satisfy g <= 0
 *   for this specialized no-slack version.
 */
[[nodiscard]] Eigen::ArrayXf mma_update(MMA_state &mma,
                                        const Eigen::ArrayXf &x,
                                        float f,
                                        const Eigen::ArrayXf &df,
                                        float g,
                                        const Eigen::ArrayXf &dg)
{
    assert(x.size() == df.size());
    assert(x.size() == dg.size());

    if (g > 1e-9f)
    {
        throw std::runtime_error(
            "mma_update requires a feasible current design");
    }

    // Convert physical design to normalized [0, 1] variables
    Eigen::ArrayXf s {(x - mma.x_min) / mma.range};
    s = s.max(0.0f).min(1.0f);

    // Gradients with respect to normalized variables
    const Eigen::ArrayXf df_ds {df * mma.range};
    const Eigen::ArrayXf dg_ds {dg * mma.range};

    // First/second iteration: use initial asymptotes
    if (mma.iteration <= 1)
    {
        mma.lower = s - mma.options.asymptote_init;
        mma.upper = s + mma.options.asymptote_init;
    }
    else
    {
        mma_update_asymptotes(mma, s);
    }

    // Subproblem bounds
    Eigen::ArrayXf alpha(mma.size);
    Eigen::ArrayXf beta(mma.size);

    for (Eigen::Index i {0}; i < mma.size; ++i)
    {
        const auto lower_bound =
            mma.lower[i] + mma.options.albefa * (s[i] - mma.lower[i]);
        const auto move_lower = s[i] - mma.options.move;

        alpha[i] = std::max({lower_bound, move_lower, 0.0f});

        const auto upper_bound =
            mma.upper[i] - mma.options.albefa * (mma.upper[i] - s[i]);
        const auto move_upper = s[i] + mma.options.move;

        beta[i] = std::min({upper_bound, move_upper, 1.0f});
    }

    // MMA reciprocal approximation
    //
    // p0, q0: objective
    // p1, q1: constraint

    Eigen::ArrayXf p0(mma.size);
    Eigen::ArrayXf q0(mma.size);
    Eigen::ArrayXf p1(mma.size);
    Eigen::ArrayXf q1(mma.size);

    for (Eigen::Index i {0}; i < mma.size; ++i)
    {
        const auto du = mma.upper[i] - s[i];
        const auto dl = s[i] - mma.lower[i];

        if (du <= 0.0f || dl <= 0.0f)
        {
            throw std::runtime_error("Invalid MMA asymptote");
        }

        const auto pos_f = std::max(df_ds[i], 0.0f);
        const auto neg_f = std::max(-df_ds[i], 0.0f);
        const auto pos_g = std::max(dg_ds[i], 0.0f);
        const auto neg_g = std::max(-dg_ds[i], 0.0f);

        // Standard MMA regularization:
        //
        // pq = 0.001 * (positive + negative) + raa0 / (xmax - xmin)
        //
        // Here normalized variable range = 1
        const auto pq_f = 0.001f * (pos_f + neg_f) + mma.options.raa0;
        const auto pq_g = 0.001f * (pos_g + neg_g) + mma.options.raa0;

        p0[i] = (pos_f + pq_f) * du * du;
        q0[i] = (neg_f + pq_f) * dl * dl;
        p1[i] = (pos_g + pq_g) * du * du;
        q1[i] = (neg_g + pq_g) * dl * dl;
    }

    // Approximation constants, ensures:
    // f_tilde(s) = f(s_current)
    // g_tilde(s) = g(s_current)

    auto r0 = f;
    auto r1 = g;

    for (Eigen::Index i {0}; i < mma.size; ++i)
    {
        const auto du = mma.upper[i] - s[i];
        const auto dl = s[i] - mma.lower[i];
        r0 -= p0[i] / du + q0[i] / dl;
        r1 -= p1[i] / du + q1[i] / dl;
    }

    // For a fixed lambda, minimize the separable Lagrangian
    const auto s_of_lambda = [&](float lambda) -> Eigen::ArrayXf
    {
        Eigen::ArrayXf result(mma.size);

        for (Eigen::Index i {0}; i < mma.size; ++i)
        {
            const auto P = p0[i] + lambda * p1[i];
            const auto Q = q0[i] + lambda * q1[i];
            const auto sqrt_P = std::sqrt(P);
            const auto sqrt_Q = std::sqrt(Q);

            // Unconstrained minimizer:
            // sqrt(P) / (U - s) = sqrt(Q) / (s - L)
            auto si = (sqrt_P * mma.lower[i] + sqrt_Q * mma.upper[i]) /
                      (sqrt_P + sqrt_Q);
            si = std::max(alpha[i], std::min(beta[i], si));

            result[i] = si;
        }

        return result;
    };

    // Evaluate approximated global constraint
    const auto g_tilde = [&](float lambda) -> float
    {
        const auto s_trial = s_of_lambda(lambda);

        auto value = r1;
        for (Eigen::Index i {0}; i < mma.size; ++i)
        {
            value += p1[i] / (mma.upper[i] - s_trial[i]);
            value += q1[i] / (s_trial[i] - mma.lower[i]);
        }

        return value;
    };

    float g_check = r1;
    for (int j = 0; j < mma.size; ++j)
    {
        g_check +=
            p1[j] / (mma.upper[j] - s[j]) + q1[j] / (s[j] - mma.lower[j]);
    }

    std::cout << "g actual      = " << g << '\n'
              << "g_tilde(s)    = " << g_check << '\n'
              << "difference    = " << g_check - g << '\n';

    float max_lower_violation = 0.0;
    float max_upper_violation = 0.0;
    for (int j = 0; j < mma.size; ++j)
    {
        max_lower_violation = std::max(max_lower_violation, alpha[j] - s[j]);

        max_upper_violation = std::max(max_upper_violation, s[j] - beta[j]);
    }

    std::cout << "max lower violation = " << max_lower_violation << '\n';
    std::cout << "max upper violation = " << max_upper_violation << '\n';

    float worst_low = 0.0;
    float worst_upp = 0.0;
    for (int j = 0; j < mma.size; ++j)
    {
        worst_low = std::max(worst_low, mma.lower[j] - alpha[j]);

        worst_upp = std::max(worst_upp, beta[j] - mma.upper[j]);
    }
    float min_gap_low = std::numeric_limits<float>::infinity();
    float min_gap_upp = std::numeric_limits<float>::infinity();
    for (int j = 0; j < mma.size; ++j)
    {
        min_gap_low = std::min(min_gap_low, s[j] - mma.lower[j]);

        min_gap_upp = std::min(min_gap_upp, mma.upper[j] - s[j]);
    }
    std::cout << "worst low = " << worst_low << '\n';
    std::cout << "worst upp = " << worst_upp << '\n';
    std::cout << "min gap low = " << min_gap_low << '\n';
    std::cout << "min gap upp = " << min_gap_upp << '\n';

    double g_terms = 0.0;

    for (int j = 0; j < mma.size; ++j)
    {
        g_terms +=
            p1[j] / (mma.upper[j] - s[j]) + q1[j] / (s[j] - mma.lower[j]);
    }

    std::cout << "|r1|       = " << std::abs(r1) << '\n'
              << "|g_terms|  = " << std::abs(g_terms) << '\n'
              << "|g|        = " << std::abs(g) << '\n';

    // First check the unconstrained solution lambda = 0
    float lambda {0.0f};

    if (g_tilde(0.0f) > 0.0f)
    {
        // Constraint active.
        // Find lambda_hi such that g_tilde(lambda_hi) <= 0

        float lambda_lo {0.0f};
        float lambda_hi {1.0f};
        bool bracketed {false};

        std::cout << "Bracketing\n";
        for (int i {0}; i < mma.options.max_expand; ++i)
        {
            std::cout << lambda_hi << " -> " << s_of_lambda(lambda_hi).mean()
                      << ' ' << s_of_lambda(lambda_hi).abs().maxCoeff() << ' '
                      << g_tilde(lambda_hi) << '\n';
            if (g_tilde(lambda_hi) <= 0.0f)
            {
                bracketed = true;
                std::cout << "Bracketed\n";
                break;
            }

            lambda_hi *= 2.0f;
        }

        if (!bracketed)
        {
            throw std::runtime_error("Could not bracket MMA dual variable");
        }

        // Scalar bisection
        for (int iter {0}; iter < mma.options.max_bisection; ++iter)
        {
            const auto lambda_mid = 0.5f * (lambda_lo + lambda_hi);
            if (g_tilde(lambda_mid) > 0.0f)
            {
                lambda_lo = lambda_mid;
            }
            else
            {
                lambda_hi = lambda_mid;
            }

            if (lambda_hi - lambda_lo <=
                mma.options.lambda_tol * (1.0f + lambda_hi))
            {
                break;
            }
        }

        lambda = 0.5f * (lambda_lo + lambda_hi);
    }

    // Recover the physical variables
    const Eigen::ArrayXf s_new {s_of_lambda(lambda)};
    const Eigen::ArrayXf x_new {mma.x_min + s_new * mma.range};

    return x_new;
}

void mma_accept(MMA_state &mma, const Eigen::ArrayXf &x)
{
    const Eigen::ArrayXf s {(x - mma.x_min) / mma.range};

    if (mma.iteration == 0)
    {
        mma.s_old_2 = s;
        mma.s_old_1 = s;
    }
    else
    {
        mma.s_old_2 = mma.s_old_1;
        mma.s_old_1 = s;
    }

    ++mma.iteration;
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
    state.areas.resize(state.elements.size());
    std::ranges::fill(state.areas, max_volume / total_length * 0.5f);

    state.sparsity_stale = true;

    const auto num_areas = static_cast<Eigen::Index>(state.areas.size());
    const auto num_coordinates = static_cast<Eigen::Index>(
        state.nodes.size() * 2 - state.immovable_dofs.size());
    const auto num_design_variables = num_areas + num_coordinates;
    Eigen::ArrayXf x_min(num_design_variables);
    Eigen::ArrayXf x_max(num_design_variables);
    x_min.block(0, 0, num_areas, 1) = min_area;
    x_max.block(0, 0, num_areas, 1) = max_area;
    x_min.block(num_areas, 0, num_coordinates, 1) = min_coordinate;
    x_max.block(num_areas, 0, num_coordinates, 1) = max_coordinate;
    mma_init(state.mma, x_min, x_max, MMA_options {});
}

void compute_activations(Optimization_state &state)
{
    // Compute activations by equal stress projection, i.e. the activations that
    // would cause the structure to be fully equally stressed.

    assert(state.element_lengths.size() == state.areas.size());

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
        const auto stiffness_constant =
            state.areas[e] * young_modulus / state.element_lengths[e];
        state.axial_forces[e] = stiffness_constant * axial_extension;
    }

    float sum {0.0f};
    for (std::size_t e {0}; e < state.areas.size(); ++e)
    {
        sum += state.element_lengths[e] * std::abs(state.axial_forces[e]);
    }
    const auto factor = max_volume / sum;

    for (std::size_t e {0}; e < state.areas.size(); ++e)
    {
        state.areas[e] = std::clamp(
            std::abs(state.axial_forces[e]) * factor, min_area, max_area);
    }
}

void geometry_step(Optimization_state &state)
{
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
        auto step = -gamma * state.dC_dx[i];
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

        optimization_init(nodes,
                          elements,
                          fixed_dofs,
                          loaded_dofs,
                          loads,
                          immovable_dofs,
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

        const auto dof =
            [&node](std::uint32_t iy, std::uint32_t ix, std::uint32_t axis)
        { return node(iy, ix) * 2 + axis; };

        const std::vector<std::uint32_t> fixed_dofs {dof(num_y / 4, 0, 0),
                                                     dof(num_y / 4, 0, 1),
                                                     dof(num_y * 3 / 4, 0, 0),
                                                     dof(num_y * 3 / 4, 0, 1)};
        const std::vector<std::uint32_t> loaded_dofs {
            dof(num_y / 2, num_x - 1, 0), dof(num_y / 2, num_x - 1, 1)};
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

        optimization_init(nodes,
                          elements,
                          fixed_dofs,
                          loaded_dofs,
                          loads,
                          immovable_dofs,
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
    /*
    Full optimization procedure:
    - Initial problem setup
    - Loop:
        - If topology changed:
            - Rebuild sparsity
            - Symbolically decompose K

        - Assemble K (only update values, same sparsity)
        - Numerically decompose K
        - Solve linear system
        - Compute sensitivities (gradients of compliance w.r.t activations and
    node positions)
        - Update design variables (MMA?)
        - If necessary, from time to time, update topology
    */

    if (state.sparsity_stale)
    {
        rebuild_stiffness_matrix_sparsity(state);
        symbolically_decompose_stiffness_matrix(state);
        state.sparsity_stale = false;
    }

    assemble_stiffness_matrix(state);
    numerically_decompose_stiffness_matrix(state);
    solve_equilibrium_system(state);
    compute_gradients(state);

    const auto num_areas = static_cast<Eigen::Index>(state.areas.size());
    const auto num_coordinates = static_cast<Eigen::Index>(
        state.nodes.size() * 2 - state.immovable_dofs.size());
    const auto num_design_variables = num_areas + num_coordinates;
    const auto C = state.compliance;
    const auto g = state.volume / max_volume - 1.0f;
    Eigen::ArrayXf x(num_design_variables);
    Eigen::ArrayXf dC(num_design_variables);
    Eigen::ArrayXf dg(num_design_variables);
    for (Eigen::Index i {0}; i < num_areas; ++i)
    {
        x[i] = state.areas[i];
        dC[i] = state.dC_dA[i];
        dg[i] = state.dV_dA[i];
    }
    std::uint32_t immovable_dof_index {0};
    std::uint32_t dof_index {0};
    for (std::uint32_t i {0}; i < state.nodes.size(); ++i)
    {
        if (immovable_dof_index < state.immovable_dofs.size() &&
            i * 2 == state.immovable_dofs[immovable_dof_index])
        {
            ++immovable_dof_index;
        }
        else
        {
            x[num_areas + dof_index] = state.nodes[i].x;
            dC[num_areas + dof_index] = state.dC_dx[i].x;
            dg[num_areas + dof_index] = state.dV_dx[i].x;
            ++dof_index;
        }

        if (immovable_dof_index < state.immovable_dofs.size() &&
            i * 2 + 1 == state.immovable_dofs[immovable_dof_index])
        {
            ++immovable_dof_index;
        }
        else
        {
            x[num_areas + dof_index] = state.nodes[i].y;
            dC[num_areas + dof_index] = state.dC_dx[i].y;
            dg[num_areas + dof_index] = state.dV_dx[i].y;
            ++dof_index;
        }
    }

    float max_element = -1e7;
    for (const auto &v : state.dC_dx)
    {
        if (const auto n = norm(v); n > max_element)
        {
            max_element = n;
        }
    }
    std::cout << "min L = " << *std::ranges::min_element(state.element_lengths)
              << ", max |dC/dx| = " << max_element << '\n';

    auto x_candidate = mma_update(state.mma, x, C, dC, g, dg);

    const auto volume_constraint = [&](const Eigen::ArrayXf &x)
    {
        std::vector<float> areas(num_areas);
        std::vector<vec2> nodes {state.nodes};

        for (Eigen::Index i {0}; i < num_areas; ++i)
        {
            areas[i] = x_candidate[i];
        }
        immovable_dof_index = 0;
        dof_index = 0;
        for (std::uint32_t i {0}; i < state.nodes.size(); ++i)
        {
            if (immovable_dof_index < state.immovable_dofs.size() &&
                i * 2 == state.immovable_dofs[immovable_dof_index])
            {
                ++immovable_dof_index;
            }
            else
            {
                nodes[i].x = x_candidate[num_areas + dof_index];
                ++dof_index;
            }

            if (immovable_dof_index < state.immovable_dofs.size() &&
                i * 2 + 1 == state.immovable_dofs[immovable_dof_index])
            {
                ++immovable_dof_index;
            }
            else
            {
                nodes[i].y = x_candidate[num_areas + dof_index];
                ++dof_index;
            }
        }

        float volume {0.0f};
        for (std::size_t e {0}; e < state.elements.size(); ++e)
        {
            const auto [node_i, node_j] = state.elements[e];
            const auto length = norm(nodes[node_j] - nodes[node_i]);
            volume += areas[e] * length;
        }

        return volume / max_volume - 1.0f;
    };

    const auto backtrack = [&]()
    {
        const auto g_candidate = volume_constraint(x_candidate);
        if (g_candidate <= 0.0f)
        {
            return x_candidate;
        }

        float alpha {1.0f};
        while (alpha > 1e-6f)
        {
            alpha *= 0.5f;

            const Eigen::ArrayXf x_trial {x + alpha * (x_candidate - x)};
            const auto g_trial = volume_constraint(x_trial);

            if (g_trial <= 0.0f)
            {
                return x_trial;
            }
        }

        return x;
    };
    const auto x_new = backtrack();

    mma_accept(state.mma, x_new);

    std::cout << std::format(
                     "[{} {}] [{} {}]    {} {}\n",
                     x_new.block(0, 0, num_areas, 1).minCoeff(),
                     x_new.block(0, 0, num_areas, 1).maxCoeff(),
                     x_new.block(num_areas, 0, num_coordinates, 1).minCoeff(),
                     x_new.block(num_areas, 0, num_coordinates, 1).maxCoeff(),
                     (x_new.block(0, 0, num_areas, 1) -
                      x.block(0, 0, num_areas, 1))
                         .abs()
                         .maxCoeff(),
                     (x_new.block(num_areas, 0, num_coordinates, 1) -
                      x.block(num_areas, 0, num_coordinates, 1))
                         .abs()
                         .maxCoeff())
              << std::flush;

    for (Eigen::Index i {0}; i < num_areas; ++i)
    {
        state.areas[i] = x_new[i];
    }
    immovable_dof_index = 0;
    dof_index = 0;
    for (std::uint32_t i {0}; i < state.nodes.size(); ++i)
    {
        if (immovable_dof_index < state.immovable_dofs.size() &&
            i * 2 == state.immovable_dofs[immovable_dof_index])
        {
            ++immovable_dof_index;
        }
        else
        {
            state.nodes[i].x = x_new[num_areas + dof_index];
            ++dof_index;
        }

        if (immovable_dof_index < state.immovable_dofs.size() &&
            i * 2 + 1 == state.immovable_dofs[immovable_dof_index])
        {
            ++immovable_dof_index;
        }
        else
        {
            state.nodes[i].y = x_new[num_areas + dof_index];
            ++dof_index;
        }
    }

    // throw std::runtime_error("");

    // compute_activations(state);
    // geometry_step(state);
}
