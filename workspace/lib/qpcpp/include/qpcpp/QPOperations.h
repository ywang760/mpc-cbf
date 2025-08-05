//
// Created by lishuo on 8/21/24.
//

#ifndef MPC_QPOPERATIONS_H
#define MPC_QPOPERATIONS_H

#include <math/Types.h>

namespace qpcpp {
template <typename T>
class QPOperations {
  public:
    using Matrix = math::Matrix<T>;
    using Vector = math::Vector<T>;
    using Row = math::Row<T>;

    class CostAddition {
      public:
        CostAddition(const Matrix& quadratic_term, const Vector& linear_term, T constant);

        const Matrix& quadratic_term() const;
        const Vector& linear_term() const;
        T constant() const;

      private:
        Matrix quadratic_term_;
        Vector linear_term_;
        T constant_;
    };

    class LinearConstraint {
      public:
        LinearConstraint(const Row& coefficients, T lower_bound, T upper_bound);

        const Row& coefficients() const;
        T lower_bound() const;
        T upper_bound() const;

      private:
        Row coefficients_;
        T lower_bound_;
        T upper_bound_;
    };

    class DecisionVariableBounds {
      public:
        DecisionVariableBounds(const Vector& lower_bounds, const Vector& upper_bounds);

        const Vector& lower_bounds() const;
        const Vector& upper_bounds() const;

      private:
        Vector lower_bounds_;
        Vector upper_bounds_;
    };
};

} // namespace qpcpp

#endif //MPC_QPOPERATIONS_H
