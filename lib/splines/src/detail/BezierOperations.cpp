//
// Created by lishuo on 2/3/24.
//

#include <math/Combinatorics.h>
#include <splines/detail/BezierOperations.h>

namespace splines {

    template <typename T>
    math::Row<T> bernsteinBasis(uint64_t bezier_degree,
                   T max_parameter, T parameter,
                   uint64_t derivative_degree) {
        if (parameter < 0 || parameter > max_parameter) {
            throw std::runtime_error("bernsteinBasis:: given parameter is outside of [0, "
                                     "max_parameter].");
        }

        math::Row<T> result = math::Row<T>::Zero(bezier_degree + 1);

        if (max_parameter == 0) {
            if (derivative_degree == 0) {
                result(0) = 1.0;
            }
            return result;
        }

        const T oneOverA = T(1.0) / max_parameter;
        for (uint64_t i = 0; i <= bezier_degree; i++) {
            T base = 0.0;
            T mult = 1.0;
            for (uint64_t j = 0; j + derivative_degree <= bezier_degree;
                 j++, mult *= parameter) {
                if (j + derivative_degree >= i) {
                    const uint64_t comb_result =
                            math::comb(bezier_degree - i, j + derivative_degree - i);

                    const uint64_t perm_result =
                            math::perm(j + derivative_degree, derivative_degree);

                    const T pow_result =
                            math::pow(oneOverA, j + derivative_degree);

                    base += (comb_result) * (pow_result) *
                            (perm_result) * mult *
                            ((j + derivative_degree - i) % 2 == 0 ? 1 : -1);
                }
            }
            const uint64_t comb_result =
                    math::comb(bezier_degree, i);

            base *= comb_result;
            result(i) = base;
        }
        return result;
    }

    /// Untested, Yutong please test this function.
    template <typename T>
    math::Matrix<T> bernsteinCoefficientMatrix(uint64_t bezier_degree, T max_parameter, uint64_t derivative_degree) {
        if (bezier_degree == std::numeric_limits<uint64_t>::max()) {
            throw std::runtime_error("bernsteinCoefficientMatrix: bezier_degree maxed, number of control points overflow.");
        }

        if (max_parameter < 0) {
            throw std::runtime_error("bernsteinCoefficientMatrix: max_parameter is negative.");
        }

        math::Matrix<T> bernstein_matrix(bezier_degree + 1, bezier_degree + 1);
        bernstein_matrix.setZero();

        if (max_parameter == 0) {
            if (derivative_degree == 0) {
                bernstein_matrix(0, 0) = 1.0;
            }
            return bernstein_matrix;
        }

        uint64_t dcombi = 1;
        const T oneOverA = T(1.0) / max_parameter;
        for (uint64_t i = 0; i < bezier_degree + 1; ++i) {
            uint64_t dminicombjmini = 1;
            T min1 = 1;

            T one_over_a_pow_j =
                    math::pow<T>(oneOverA, i);

            for (uint64_t j = i; j < bezier_degree + 1;
                 ++j, min1 *= -1, one_over_a_pow_j *= oneOverA) {
                bernstein_matrix(i, j) =
                        dcombi * dminicombjmini * min1 * one_over_a_pow_j;

                if (bezier_degree != j &&
                    dminicombjmini > std::numeric_limits<uint64_t>::max() /
                                     (bezier_degree - j)) {
                    throw std::runtime_error("bernsteinCoefficientMatrix: dminicombjmini overflows.");
                }

                dminicombjmini *= (bezier_degree - j);
                dminicombjmini /= (j + 1 - i);
            }

            if ((bezier_degree != i) &&
                dcombi >
                std::numeric_limits<uint64_t>::max() / (bezier_degree - i)) {
                throw std::runtime_error("bernsteinCoefficientMatrix: dcombi overflows.");
            }

            dcombi *= (bezier_degree - i);
            dcombi /= (i + 1);
        }

        math::Matrix<T> derivative(bezier_degree + 1, bezier_degree + 1);
        derivative.setZero();

        uint64_t jpermk =
                math::fac(derivative_degree);

        for (uint64_t j = derivative_degree; j < bezier_degree + 1; ++j) {
            derivative(j, j - derivative_degree) = jpermk;

            if (jpermk > std::numeric_limits<uint64_t>::max() / (j + 1)) {
                throw std::runtime_error("bernsteinCoefficientMatrix: jpermk overflows.");
            }

            jpermk *= (j + 1);
            jpermk /= (j + 1 - derivative_degree);
        }

        return bernstein_matrix * derivative;
    }

    template math::Row<double> bernsteinBasis(uint64_t bezier_degree,
            double max_parameter, double parameter,
            uint64_t derivative_degree);

    template math::Row<float> bernsteinBasis(uint64_t bezier_degree,
            float max_parameter, float parameter,
            uint64_t derivative_degree);

    template math::Matrix<double> bernsteinCoefficientMatrix(uint64_t bezier_degree,
            double max_parameter,
            uint64_t derivative_degree);

    template math::Matrix<float> bernsteinCoefficientMatrix(uint64_t bezier_degree,
            float max_parameter,
            uint64_t derivative_degree);

} // splines