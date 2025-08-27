#include "analysis.hpp"

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

    for (const auto [i, j] : state.elements)
    {
        const auto vec = state.nodes[j] - state.nodes[i];
        const auto length = norm(vec);
        const auto [c, s] = vec / length;
        const auto EA_over_length = (young_modulus * area) / length;

        const float element_stiffness_matrix[4][4] {
            {c * c, c * s, -c * c, -c * s},
            {c * s, s * s, -c * s, -s * s},
            {-c * c, -c * s, c * c, c * s},
            {-c * s, -s * s, c * s, s * s}};

        const int dof_indices[4] {dof_all_to_free[2 * i],
                                  dof_all_to_free[2 * i + 1],
                                  dof_all_to_free[2 * j],
                                  dof_all_to_free[2 * j + 1]};
        for (unsigned int a {0}; a < 4; ++a)
        {
            for (unsigned int b {0}; b < 4; ++b)
            {
                if (dof_indices[a] != -1 && dof_indices[b] != -1)
                {
                    triplets.emplace_back(dof_indices[a],
                                          dof_indices[b],
                                          EA_over_length *
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

    state.displacements.setZero(num_dofs);
}

void solve(Analysis_state &state)
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
    state.displacements.setZero();
    state.displacements(state.free_dofs) = free_displacements;
}
