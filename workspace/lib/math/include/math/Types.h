#ifndef MATH_TYPES_H
#define MATH_TYPES_H

#include <Eigen/Dense>

namespace math {

template <typename T, unsigned int DIM>
using AlignedBox = Eigen::AlignedBox<T, DIM>;

template <typename T, unsigned int DIM>
using Hyperplane = Eigen::Hyperplane<T, DIM>;

template <typename T, unsigned int DIM>
using VectorDIM = Eigen::Matrix<T, DIM, 1>;

template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T, unsigned int R, unsigned int C>
using MatrixRC = Eigen::Matrix<T, R, C>;

template <typename T, unsigned int DIM>
using CholeskyDecomposition = Eigen::LLT<MatrixRC<T, DIM, DIM>>;
} // namespace math

#endif //MATH_TYPES_H
