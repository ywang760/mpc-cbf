//
// Created by lishuo on 8/31/24.
//

#include <cbf/optimization/CBFQPOperations.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    CBFQPOperations<T, DIM>::CBFQPOperations(std::shared_ptr<FovCBF> cbf) : cbf_(cbf) {}

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::CostAddition
    CBFQPOperations<T, DIM>::desiredControlCost(const VectorDIM &desired_u) {
        Matrix quadratic_term(DIM, DIM);
        quadratic_term.setZero();
        Vector linear_term(DIM);
        linear_term.setZero();

        // set the quadratic term
        quadratic_term.setIdentity();
        // set the linear term
        linear_term = -2.0 * desired_u;
        return CostAddition(quadratic_term, linear_term, 0);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::safetyConstraint(const Vector &state,
                                              const Vector &target_state) {
        Vector coefficients = cbf_->getSafetyConstraints(state, target_state);
        T bound = cbf_->getSafetyBound(state, target_state);
        return LinearConstraint(coefficients, -bound, std::numeric_limits<T>::max());
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::leftBorderConstraint(const Vector &state,
                                                  const Vector &target_state) {
        Vector coefficients = cbf_->getLBConstraints(state, target_state);
        T bound = cbf_->getLBBound(state, target_state);
        return LinearConstraint(coefficients, -bound, std::numeric_limits<T>::max());
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rightBorderConstraint(const Vector &state,
                                                   const Vector &target_state) {
        Vector coefficients = cbf_->getRBConstraints(state, target_state);
        T bound = cbf_->getRBBound(state, target_state);
        return LinearConstraint(coefficients, -bound, std::numeric_limits<T>::max());
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rangeConstraint(const Vector &state,
                                             const Vector &target_state) {
        Vector coefficients = cbf_->getRangeConstraints(state, target_state);
        T bound = cbf_->getRangeBound(state, target_state);
        return LinearConstraint(coefficients, -bound, std::numeric_limits<T>::max());
    }

    template class CBFQPOperations<double, 3U>;

} // cbf