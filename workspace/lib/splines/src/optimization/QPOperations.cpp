//
// Created by lishuo on 8/21/24.
//

#include <splines/optimization/QPOperations.h>

namespace splines {
    /// Define the objective function form,
    /// it contains quadratic, linear, constant terms.
    template <typename T>
    QPOperations<T>::CostAddition::CostAddition(
            const Matrix &quadratic_term, const Vector &linear_term, T constant)
            : quadratic_term_(quadratic_term),
              linear_term_(linear_term),
              constant_(constant) {}

    template <typename T>
    const typename QPOperations<T>::Matrix &
    QPOperations<T>::CostAddition::quadratic_term() const {
        return quadratic_term_;
    }

    template <typename T>
    const typename QPOperations<T>::Vector &
    QPOperations<T>::CostAddition::linear_term() const {
        return linear_term_;
    }

    template <typename T>
    T QPOperations<T>::CostAddition::constant() const {
        return constant_;
    }

    /// Define the linear constraint form,
    /// it has the coefficients, lower_bound, upper_bound terms.
    template <typename T>
    QPOperations<T>::LinearConstraint::LinearConstraint(
            const Row &coefficients, T lower_bound, T upper_bound)
            : coefficients_(coefficients),
              lower_bound_(lower_bound),
              upper_bound_(upper_bound) {}

    template <typename T>
    const typename QPOperations<T>::Row &
    QPOperations<T>::LinearConstraint::coefficients() const {
        return coefficients_;
    }

    template <typename T>
    T QPOperations<T>::LinearConstraint::lower_bound() const {
        return lower_bound_;
    }

    template <typename T>
    T QPOperations<T>::LinearConstraint::upper_bound() const {
        return upper_bound_;
    }

    /// Define the decisionVariableBounds form,
    /// it has lower_bounds and upper_bounds.
    template <typename T>
    QPOperations<T>::DecisionVariableBounds::DecisionVariableBounds(
            const Vector &lower_bounds, const Vector &upper_bounds)
            : lower_bounds_(lower_bounds),
              upper_bounds_(upper_bounds) {}

    template <typename T>
    const typename QPOperations<T>::Vector &
    QPOperations<T>::DecisionVariableBounds::lower_bounds() const {
        return lower_bounds_;
    }

    template <typename T>
    const typename QPOperations<T>::Vector &
    QPOperations<T>::DecisionVariableBounds::upper_bounds() const {
        return upper_bounds_;
    }

    template class QPOperations<float>::CostAddition;
    template class QPOperations<double>::CostAddition;

    template class QPOperations<float>::LinearConstraint;
    template class QPOperations<double>::LinearConstraint;

    template class QPOperations<float>::DecisionVariableBounds;
    template class QPOperations<double>::DecisionVariableBounds;
} // splines