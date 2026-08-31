#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "vec.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <array>
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
    std::vector<std::uint32_t> free_dofs;
    std::vector<std::uint32_t> all_to_free_dofs;
    Eigen::SparseMatrix<float, Eigen::ColMajor> stiffness_matrix;
    std::vector<std::array<std::uint32_t, 16>> element_k_indices;
    std::vector<vec2> element_directions;
    std::vector<float> element_lengths;
    std::vector<float> activations;
    Eigen::SimplicialLDLT<decltype(stiffness_matrix), Eigen::Lower> solver;
    Eigen::VectorXf loads;
    Eigen::VectorXf displacements;
    bool sparsity_stale;

    std::vector<std::uint32_t> movable_nodes;
    std::vector<float> new_activations;
    std::vector<float> axial_forces;
    std::vector<vec2> dC_dx;
    float compliance;
    float predicted_compliance_reduction;
    float delta_x;
    int step;
};

enum struct Problem
{
    regular_grid,
    random_delaunay,
    hexagonal
};

void optimization_create_problem(Optimization_state &state, Problem problem);
void optimization_step(Optimization_state &state);

#endif