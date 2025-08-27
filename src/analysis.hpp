#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "vec.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstdint>
#include <vector>

struct Element
{
    std::uint32_t node_i;
    std::uint32_t node_j;
};

struct Analysis_state
{
    std::vector<vec2> nodes;
    std::vector<Element> elements;
    std::vector<float> activations;
    std::vector<float> stiffness_constants;
    std::vector<vec2> element_directions;
    std::vector<std::uint32_t> fixed_dofs;
    std::vector<std::uint32_t> free_dofs;
    Eigen::SparseMatrix<float> stiffness_matrix;
    Eigen::VectorXf loads;
};

struct Linear_system_result
{
    Eigen::VectorXf displacements;
    std::vector<float> axial_forces;
    std::vector<float> energies;
};

void assemble(Analysis_state &state);
[[nodiscard]] Linear_system_result solve(const Analysis_state &state);

void make_triangulation(Analysis_state &state);

#endif