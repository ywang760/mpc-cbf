//
// Created by lishuo on 2/13/24.
//

#include <math/Combinatorics.h>
#include <splines/curves/Bezier.h>
#include <splines/optimization/BezierQPOperations.h>

namespace splines {
    template <typename T, unsigned int DIM>
    BezierQPOperations<T, DIM>::BezierQPOperations(Params &p)
    : num_control_points_(p.num_control_points_),
    max_parameter_(p.max_parameter_) {}

    template <typename T, unsigned int DIM>
    std::size_t BezierQPOperations<T, DIM>::numDecisionVariables() const {
        return DIM * num_control_points_;
    }

    template <typename T, unsigned int DIM>
    std::unique_ptr<typename BezierQPOperations<T, DIM>::SingleParameterCurve>
    BezierQPOperations<T, DIM>::generateCurveFromSolution(const Vector &decision_variables) const {
        using Bezier = splines::Bezier<T, DIM>;
        // pre-condition throw
        if (decision_variables.rows() != numDecisionVariables()) {
            throw std::runtime_error("number of decision variables does not match the required number.");
        }
        // create bezier curve
        std::unique_ptr<Bezier> bezier_ptr = std::make_unique<Bezier>(max_parameter_, std::vector<VectorDIM>());

        // append the control points from decision variable
        for (std::size_t control_point_idx = 0; control_point_idx < num_control_points_; ++control_point_idx) {
            VectorDIM control_point;
            for (unsigned int d = 0; d < DIM; ++d) {
                control_point(d) = decision_variables(d * num_control_points_ + control_point_idx);
            }
            bezier_ptr->appendControlPoint(control_point);
        }

        return std::unique_ptr<SingleParameterCurve>(bezier_ptr.release());
    }

    /// example of decision variables:
    /// for 5 control points:
    /// 0,1,2,3,4, 5,6,7,8,9, 10,11,12,13,14
    /// -> 0,5,10 .. 1,6,11 .. 2,7,12 .. 3,8,13 .. 4,9,14
    template <typename T, unsigned int DIM>
    typename BezierQPOperations<T, DIM>::Row
    BezierQPOperations<T, DIM>::evalBasisRow(unsigned int dimension, T parameter, uint64_t derivative_degree) const {
        if (dimension >= DIM) {
            throw std::runtime_error("dimension is out of range of DIM");
        }

        if (num_control_points_ == 0) {
            Row result(0);
            return result;
        }
        // compute the bernstein basis
        Row basis = bernsteinBasis(num_control_points_ - 1, max_parameter_, parameter, derivative_degree);

        Row result = Row::Zero(numDecisionVariables());
        result.block(/*start_row=*/0, /*start_col=*/dimension * num_control_points_,
                     /*num_rows=*/1, /*num_cols=*/num_control_points_) = basis;
        return result;
    }

    template <typename T, unsigned int DIM>
    typename BezierQPOperations<T, DIM>::CostAddition
    BezierQPOperations<T, DIM>::integratedSquaredDerivativeCost(uint64_t derivative_degree, T lambda) const {
//        using Matrix = BezierQPOperations<T, DIM>::Matrix;
        Matrix quadratic_term(numDecisionVariables(), numDecisionVariables());
        quadratic_term.setZero();
        Vector linear_term(numDecisionVariables());
        linear_term.setZero();

        // if no control points, return 0 cost.
        if (num_control_points_ == 0) {
            return CostAddition(quadratic_term, linear_term, 0);
        }

        if (derivative_degree <= /*degree=*/num_control_points_ - 1) {
            const Matrix bernstein_coefficient_matrix =
                    bernsteinCoefficientMatrix(/*bezier_degree=*/num_control_points_ - 1,
                                                                 max_parameter_,
                                                                 derivative_degree);

            Matrix SQI(num_control_points_, num_control_points_);
            SQI.setZero();

            for (std::size_t i = 0; i < num_control_points_; i++) {
                for (std::size_t j = 0; j < num_control_points_; j++) {
                    const T pow =
                            math::pow<T>(max_parameter_, i + j + 1);

                    SQI(i, j) = pow / (i + j + 1);
                }
            }

            Matrix cost = lambda * bernstein_coefficient_matrix * SQI *
                           bernstein_coefficient_matrix.transpose();

            for (unsigned int dimension = 0; dimension < DIM; ++dimension) {
                quadratic_term.block(dimension * num_control_points_,
                                     dimension * num_control_points_,
                                     num_control_points_, num_control_points_) =
                        cost;
            }
        }

        return CostAddition(quadratic_term, linear_term, 0);
    }

    template <typename T, unsigned int DIM>
    typename BezierQPOperations<T, DIM>::CostAddition
    BezierQPOperations<T, DIM>::evalCost(T parameter, uint64_t derivative_degree,
                                         const VectorDIM &target, T lambda) const {
        Matrix quadratic_term(numDecisionVariables(), numDecisionVariables());
        quadratic_term.setZero();
        Vector linear_term(numDecisionVariables());
        linear_term.setZero();

        if (num_control_points_ == 0) {
            return CostAddition(quadratic_term, linear_term, /*constant=*/0);
        }

        const Row& basis = bernsteinBasis(num_control_points_ - 1, max_parameter_,
                                           parameter, derivative_degree);

        const Matrix quadratic_term_per_dimension =
                lambda * basis.transpose() * basis;
        const Vector linear_term_per_dimension_except_target =
                -2 * lambda * basis.transpose();
        T constant_except_lambda = 0;

        for (unsigned int dimension = 0; dimension < DIM; ++dimension) {
            quadratic_term.block(dimension * num_control_points_, dimension * num_control_points_,
                                 num_control_points_, num_control_points_) =
                    quadratic_term_per_dimension;

            linear_term.block(dimension * num_control_points_, 0, num_control_points_, 1) =
                    linear_term_per_dimension_except_target * target(dimension);

            constant_except_lambda += target(dimension) * target(dimension);
        }

        return CostAddition(quadratic_term, linear_term,
                            constant_except_lambda * lambda);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename BezierQPOperations<T, DIM>::LinearConstraint>
    BezierQPOperations<T, DIM>::evalConstraint(T parameter, uint64_t derivative_degree,
                                               const VectorDIM& target) const {
        std::vector<LinearConstraint> linear_constraints;
        if (num_control_points_ == 0) {
            return linear_constraints;
        }
        const Row &basis = bernsteinBasis(/*bezier_degree=*/num_control_points_ - 1,
        max_parameter_, parameter, derivative_degree);
        // evalConstraint for each dimension
        for (unsigned int dimension = 0; dimension < DIM; dimension++) {
            Row coefficients(numDecisionVariables());
            coefficients.setZero();
            coefficients.block(0, dimension*num_control_points_, 1, num_control_points_) = basis;
            linear_constraints.push_back(LinearConstraint(coefficients, target(dimension), target(dimension)));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    typename BezierQPOperations<T, DIM>::DecisionVariableBounds
    BezierQPOperations<T, DIM>::boundingBoxConstraint(
            const AlignedBox& bounding_box) const {

        Vector lower_bounds(numDecisionVariables());
        Vector upper_bounds(numDecisionVariables());

        for (unsigned int dimension_idx = 0; dimension_idx < DIM; ++dimension_idx) {
            T bounding_box_min_for_dimension = bounding_box.min()(dimension_idx);
            T bounding_box_max_for_dimension = bounding_box.max()(dimension_idx);
            for (std::size_t control_point_idx = 0;
                 control_point_idx < num_control_points_; ++control_point_idx) {
                lower_bounds(dimension_idx * num_control_points_ +
                             control_point_idx) = bounding_box_min_for_dimension;
                upper_bounds(dimension_idx * num_control_points_ +
                             control_point_idx) = bounding_box_max_for_dimension;
            }
        }
        return DecisionVariableBounds(lower_bounds, upper_bounds);
    }

    template <typename T, unsigned int DIM>
    typename BezierQPOperations<T, DIM>::DecisionVariableBounds
    BezierQPOperations<T, DIM>::positionConstraint(const VectorDIM &position) const {

        Vector lower_bounds(numDecisionVariables());
        Vector upper_bounds(numDecisionVariables());

        for (unsigned int dimension_idx = 0; dimension_idx < DIM; ++dimension_idx) {
            T position_for_dimension = position(dimension_idx);
            for (std::size_t control_point_idx = 0;
                 control_point_idx < num_control_points_; ++control_point_idx) {
                lower_bounds(dimension_idx * num_control_points_ +
                             control_point_idx) = position_for_dimension;
                upper_bounds(dimension_idx * num_control_points_ +
                             control_point_idx) = position_for_dimension;
            }
        }
        return DecisionVariableBounds(lower_bounds, upper_bounds);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename BezierQPOperations<T, DIM>::LinearConstraint>
    BezierQPOperations<T, DIM>::boundingBoxConstraintAll(
            const AlignedBox& bounding_box, uint64_t derivative_degree) const {
        if (math::isApproximatelyEqual<T>(max_parameter_, T(0.0))) {
            throw std::runtime_error("max_parameter_ is zero");
        }

        const VectorDIM& min = bounding_box.min();
        const VectorDIM& max = bounding_box.max();

        std::vector<LinearConstraint> constraints;

        const T perm =
                math::perm(num_control_points_ - 1, derivative_degree);

        const T pow = math::pow(T(1.0) / max_parameter_, derivative_degree);

        for (uint64_t i = 0; i < num_control_points_ - derivative_degree; ++i) {
            Row control_point_multipliers = Row::Zero(num_control_points_);
            int minusOneOverJ = 1;
            for (uint64_t j = 0; j <= derivative_degree; ++j) {
                const T comb = math::comb(derivative_degree, j);
                control_point_multipliers(i + derivative_degree - j) =
                        pow * perm * comb * minusOneOverJ;
                minusOneOverJ *= -1;
            }

            for (unsigned int d = 0; d < DIM; ++d) {
                Row coefficients = Row::Zero(numDecisionVariables());
                coefficients.block(0, d * num_control_points_, 1,
                                   num_control_points_) = control_point_multipliers;
                constraints.push_back(LinearConstraint(
                        coefficients, bounding_box.min()(d), bounding_box.max()(d)));
            }
        }

        return constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename BezierQPOperations<T, DIM>::LinearConstraint>
    BezierQPOperations<T, DIM>::hyperplaneConstraintAll(
            const Hyperplane& hyperplane, T epsilon) const {

        std::vector<LinearConstraint> linear_constraints;
        for (std::size_t control_point_idx = 0;
             control_point_idx < num_control_points_; ++control_point_idx) {
            Row coefficients(numDecisionVariables());
            coefficients.setZero();
            for (unsigned int dimension_idx = 0; dimension_idx < DIM;
                 ++dimension_idx) {
                coefficients(dimension_idx * num_control_points_ +
                             control_point_idx) =
                        hyperplane.normal()(dimension_idx);
            }
            linear_constraints.push_back(
                    LinearConstraint(coefficients, std::numeric_limits<T>::lowest(),
                                     -hyperplane.offset() - epsilon));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename BezierQPOperations<T, DIM>::LinearConstraint>
    BezierQPOperations<T, DIM>::hyperplaneConstraintAt(
            T parameter, const Hyperplane& hyperplane, T epsilon) const {

        std::vector<LinearConstraint> linear_constraints;

        if (num_control_points_ == 0) {
            return linear_constraints;
        }

        Row basis = bernsteinBasis(
                /*degree=*/num_control_points_ - 1, max_parameter_, parameter,
                /*derivative_degree=*/0);

        Row coefficients(numDecisionVariables());
        coefficients.setZero();
        for (unsigned int dimension_idx = 0; dimension_idx < DIM; ++dimension_idx) {
            coefficients.block(0, dimension_idx * num_control_points_, 1,
                               num_control_points_) =
                    basis * hyperplane.normal()(dimension_idx);
        }
        linear_constraints.push_back(
                LinearConstraint(coefficients, std::numeric_limits<T>::lowest(),
                                 -hyperplane.offset()  - epsilon));
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    void BezierQPOperations<T, DIM>::set_max_parameter(T new_max_parameter) {
        if (new_max_parameter < 0) {
            throw std::invalid_argument("BezierQPOperations::set_max_parameter: "
                                        "new_max_parameter is negative. new_max_parameter: " +
                                        std::to_string(new_max_parameter));
        }

        max_parameter_ = new_max_parameter;
    }

    template <typename T, unsigned int DIM>
    T BezierQPOperations<T, DIM>::max_parameter() const {
        return max_parameter_;
    }

    template class BezierQPOperations<float, 2U>;
    template class BezierQPOperations<double, 2U>;
    template class BezierQPOperations<float, 3U>;
    template class BezierQPOperations<double, 3U>;

} // splines