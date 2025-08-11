//
// Created on 8/4/25.
//

#ifndef MPC_CBF_MPCCBFQPGENERATORBASE_H
#define MPC_CBF_MPCCBFQPGENERATORBASE_H

#include <math/Helpers.h>
#include <memory>
#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>
#include <mpc_cbf/optimization/MPCCBFQPOperationsBase.h>
#include <vector>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class MPCCBFQPGeneratorBase {
  public:
    using Problem = qpcpp::Problem<T>;
    using State = model::State<T, DIM>;
    using Vector = math::Vector<T>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using Row = math::Row<T>;
    using LinearConstraint = typename qpcpp::QPOperations<T>::LinearConstraint;
    using CostAddition = typename qpcpp::QPOperations<T>::CostAddition;
    using PiecewiseBezierMPCQPGenerator = mpc::PiecewiseBezierMPCQPGenerator<T, DIM>;

    virtual ~MPCCBFQPGeneratorBase() = default;

    // Common interface
    Problem& problem();
    void addSlackCost(const std::vector<double>& slack_weights,
                      const std::vector<qpcpp::Variable<T>*>& specific_slack_variables = {});
    std::shared_ptr<PiecewiseBezierMPCQPGenerator> piecewise_mpc_qp_generator_ptr();

    // Common slack variable management
    void addSlackVariables(int num_neighbors);
    void addCostAdditionForSlackVariables(const CostAddition& cost_addition,
                                          const std::vector<qpcpp::Variable<T>*>& specific_slack_variables = {});
    void
    addLinearConstraintForPiecewiseWithSlackVariables(const LinearConstraint& linear_constraint,
                                                      const Row& slack_coefficients,
                                                      const std::vector<qpcpp::Variable<T>*>& specific_slack_variables = {});

  protected:
    // Constructor for derived classes
    MPCCBFQPGeneratorBase();

    // Common members
    std::shared_ptr<PiecewiseBezierMPCQPGenerator> piecewise_mpc_qp_generator_ptr_;
    std::vector<qpcpp::Variable<T>*> slack_variables_;
    bool slack_mode_;
    int num_neighbors_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_MPCCBFQPGENERATORBASE_H