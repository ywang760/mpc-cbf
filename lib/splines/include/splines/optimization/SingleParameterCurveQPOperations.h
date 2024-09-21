//
// Created by lishuo on 2/13/24.
//

#ifndef SPLINES_SINGLEPARAMETERCURVEQPOPERATIONS_H
#define SPLINES_SINGLEPARAMETERCURVEQPOPERATIONS_H

#include <qpcpp/QPOperations.h>
#include <splines/curves/SingleParameterCurve.h>

namespace splines {

    template <typename T, unsigned int DIM>
    class BezierQPOperations;

    template <typename T, unsigned int DIM>
    class SingleParameterCurveQPOperations {
    public:
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Row = math::Row<T>;
        using AlignedBox = math::AlignedBox<T, DIM>;
        using Hyperplane = math::Hyperplane<T, DIM>;
        using SingleParameterCurve = splines::SingleParameterCurve<T, DIM>;
        using CostAddition = typename qpcpp::QPOperations<T>::CostAddition;
        using LinearConstraint = typename qpcpp::QPOperations<T>::LinearConstraint;
        using DecisionVariableBounds = typename qpcpp::QPOperations<T>::DecisionVariableBounds;

        virtual ~SingleParameterCurveQPOperations() = default;

        // ----------------
        // virtual functions
        virtual std::unique_ptr<SingleParameterCurve> generateCurveFromSolution(const Vector& decision_variables) const = 0;
        virtual Row evalBasisRow(unsigned int dimension, T parameter, uint64_t derivative_degree) const = 0;

        /// objective functions
        // returns Q, q and c that adds the cost lambda * \int_{0}^{max_parameter}
        // ||df^{derivative_degree}}(u)/du^{derivative_degree}||_2^2 du as p^T Q p +
        // p^Tq + c where p denotes the decision variables. return value is not ok
        // if anything overflows during computation
        virtual CostAddition integratedSquaredDerivativeCost(uint64_t derivative_degree, T lambda) const = 0;

        // returns Q, q and c that adds the cost lambda *
        // ||f^{derivative_degree}(parameter)-target||^2 as p^T Q p + p^Tq + c where
        // p denotes the decision variables. return status is not ok if parameter is
        // our of range
        virtual CostAddition evalCost(T parameter, uint64_t derivative_degree,
                                      const VectorDIM& target, T lambda) const = 0;

        //  returns set of constraints that enforces
        //  f^{derivative_degree}(parameter) = target. return status is not ok if
        //  parameter is out of range or anything overflows.
        virtual std::vector<LinearConstraint> evalConstraint(
                T parameter, uint64_t derivative_degree,
                const VectorDIM& target) const = 0;

        // returns the number of decision variables of the QP that will be generated
        // from these operations.
        virtual std::size_t numDecisionVariables() const = 0;

        // sets the max parameter of the curve to the given value. curve is assumed
        // to be defined in [0, new_max_parameter]
        virtual void set_max_parameter(T new_max_parameter) = 0;

        // get the max parameter of the curve
        virtual T max_parameter() const = 0;

        // returns lower and upper bounds on decision variables that enforces f(u)
        // is in the given bounding box for all u
        virtual DecisionVariableBounds boundingBoxConstraint(
                const AlignedBox& bounding_box) const = 0;

        // return lower and upper bounds on the decision variables that enforces f(u)
        // is on the given position for all u
        virtual DecisionVariableBounds positionConstraint(const VectorDIM& position) const = 0;

        /**
         * @brief Constaint the given derivative of the curve to be inside the given
         * bounding box at all points
         *
         * @param bounding_box Bounding box that the given derivative of the curve
         * must be contained in
         * @param derivative_degree Derivative degree of the curve
         * @return std::vector<LinearConstraint> Set of constraints that enforces
         * the given derivative of the curve to be inside the given bounding box.
         */
        virtual std::vector<LinearConstraint>
        boundingBoxConstraintAll(const AlignedBox& bounding_box,
                                 uint64_t derivative_degree) const = 0;

        // given hyperplane hp(x) = a_1x_1+a_2x_2+...+a_nx_n + d = 0 such that
        // n=(a_1, ..., a_n) is the normal returns set of constraints that enforces
        // hp(f(u)) <= 0 for all u, add epsilon to make hp(f(u)) < 0
        virtual std::vector<LinearConstraint> hyperplaneConstraintAll(
                const Hyperplane& hyperplane, T epsilon) const = 0;

        // given hyperplane hp(x) = a_1x_1+a_2x_2+...+a_nx_n + d = 0
        // such that n=(a_1, ..., a_n) is the normal,
        // returns set of constraints that enforces hp(f(parameter)) <= 0
        virtual std::vector<LinearConstraint>
        hyperplaneConstraintAt(T parameter, const Hyperplane& hyperplane, T epsilon) const = 0;

    };

} // splines

#endif //SPLINES_SINGLEPARAMETERCURVEQPOPERATIONS_H
