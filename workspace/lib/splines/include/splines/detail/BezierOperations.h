//
// Created by lishuo on 2/3/24.
//

#ifndef SPLINES_BEZIEROPERATIONS_H
#define SPLINES_BEZIEROPERATIONS_H

#include <math/Types.h>

namespace splines {
/**
     * @brief For a Bezier curve with degree bezier_degree and max parameter
     * max_parameter, return the bernstein basis polynomials that would result in
     * the evaluation of the derivative_degree^{th} derivative of bezier curve at
     * parameter when multiplied with control points
     *
     * @tparam T Floating point number type
     * @param bezier_degree Degree of the Bezier curve
     * @param max_parameter Maximum parameter of the Bezier curve
     * @param parameter Parameter at which the derivative is to be evaluated
     * @param derivative_degree Degree of the derivative to be evaluated
     * @return absl::StatusOr<math::Row<T>> Row vector containing the evaluations of
     * bernstein basis polynomials at parameter. If parameter is outside the range
     * [0, max_parameter], return status is not OK. If anything overflows during the
     * computation, return status is not OK.
     */
template <typename T>
math::Row<T> bernsteinBasis(uint64_t bezier_degree, T max_parameter, T parameter,
                            uint64_t derivative_degree);

/**
     * @brief Compute the polynomial coefficient matrix for the
     * {derivative_degree}^{th} derivative of bernstein base functions where each
     * row r with index i contains polynomial coefficients for the
     * {derivative_degree}^{th} derivative of i^th bernstein polynomial of degree
     * bezier_degree such that the bernstein polynomial is r(0) + r(1) t + r(2) t^2
     * + ... + r(d) t^d.
     *
     * @tparam T Floating point number type
     * @param bezier_degree Degree of the Bezier curve
     * @param max_parameter Maximum parameter of the Bezier curve, and hence the
     * berstein polynomials
     * @param derivative_degree Degree of the derivative for which the coefficients
     * are evaluated
     * @return absl::StatusOr<math::Matrix<T>> Matrix containing the polynomial
     * coefficients for the {derivative_degree}^{th} derivative of bernstein base
     * functions for the Bezier curve with degree bezier_degree and max parameter =
     * max_parameter. Return status is not OK if anything overflows during the
     * computation, or bezier_degree = maximum uint64_t or max_parameter is
     * negative.
     */
template <typename T>
math::Matrix<T> bernsteinCoefficientMatrix(uint64_t bezier_degree, T max_parameter,
                                           uint64_t derivative_degree);

} // namespace splines

#endif //SPLINES_BEZIEROPERATIONS_H
