#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "vec.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <cstdint>
#include <vector>

struct Element
{
    std::uint32_t node_i;
    std::uint32_t node_j;
};

struct MMA_options
{
    // Maximum fraction of the physical variable range allowed in one MMA
    // iteration
    float move = 0.1f;

    // Initial asymptote distance, as a fraction of range
    float asymptote_init = 0.3f;

    // Asymptote adaptation
    float asymptote_inc = 1.1f;
    float asymptote_dec = 0.7f;

    // Distance between current point and asymptotes used for alpha/beta
    float albefa = 0.1f;

    // Standard MMA regularization parameter
    float raa0 = 1e-5f;

    // Bisection stopping tolerance
    float lambda_tol = 1e-9f;

    // Maximum number of lambda bracket expansions
    int max_expand = 100;

    // Maximum number of bisection iterations
    int max_bisection = 100;
};

struct MMA_state
{
    Eigen::Index size;
    Eigen::ArrayXf x_min;
    Eigen::ArrayXf x_max;
    Eigen::ArrayXf range;
    Eigen::ArrayXf s_old_1;
    Eigen::ArrayXf s_old_2;
    Eigen::ArrayXf lower;
    Eigen::ArrayXf upper;
    int iteration;
    bool initialized;
    MMA_options options;
};

struct Optimization_state
{
    std::vector<vec2> nodes;
    std::vector<Element> elements;

    std::vector<std::uint32_t> fixed_dofs;
    std::vector<std::uint32_t> free_dofs;
    std::vector<std::uint32_t> all_to_free_dofs;
    std::vector<std::uint32_t> immovable_dofs;

    std::vector<float> areas;

    // For each element, for each coefficient of the local 4x4 matrix, index
    // into stiffness_matrix.valuePtr()
    std::vector<std::array<std::uint32_t, 16>> element_k_indices;
    std::vector<vec2> element_directions;
    std::vector<float> element_lengths;
    Eigen::SparseMatrix<float, Eigen::ColMajor> stiffness_matrix;
    bool sparsity_stale;
    Eigen::SimplicialLDLT<decltype(stiffness_matrix), Eigen::Lower> solver;
    Eigen::VectorXf loads;
    Eigen::VectorXf displacements;

    std::vector<float> axial_forces;

    float compliance;
    float volume;
    std::vector<float> dC_dA;
    std::vector<vec2> dC_dx;
    std::vector<float> dV_dA;
    std::vector<vec2> dV_dx;

    MMA_state mma;
};

enum struct Problem
{
    regular_grid,
    random_delaunay
};

void optimization_create_problem(Optimization_state &state, Problem problem);
void optimization_step(Optimization_state &state);

#endif