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
    std::vector<std::uint32_t> fixed_dofs;
    std::vector<std::uint32_t> free_dofs;
    Eigen::SparseMatrix<float> stiffness_matrix;
    Eigen::VectorXf loads;
    Eigen::VectorXf displacements;
};

void assemble(Analysis_state &state);
void solve(Analysis_state &state);

#endif