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

struct Optimization_state
{
    std::vector<vec2> nodes;
    std::vector<Element> elements;
    std::vector<std::uint32_t> fixed_dofs;
    std::vector<std::uint32_t> free_dofs;
    std::vector<std::uint32_t> all_to_free_dofs;
    std::vector<std::uint32_t> immovable_dofs;

    std::vector<float> activations;

    std::vector<vec2> element_directions;
    std::vector<float> element_lengths;

    Eigen::SparseMatrix<float> stiffness_matrix;
    Eigen::VectorXf loads;
    Eigen::VectorXf displacements;

    std::vector<float> axial_forces;
    std::vector<float> energies;
    std::vector<vec2> gradients;
};

void optimization_create_problem(Optimization_state &state);
void optimization_step(Optimization_state &state);

#endif