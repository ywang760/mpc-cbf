#ifndef MPC_CONNECTIVITYQPGENERATOR_H
#define MPC_CONNECTIVITYQPGENERATOR_H

#include <cbf/detail/ConnectivityCBF.h>
#include <cbf/optimization/CBFQPGeneratorBase.h>

namespace cbf {
/**
     * @brief Connectivity Quadratic Program Generator
     * 
     * This class implements the connectivity-specific constraint methods for CBF QP generation.
     * It handles connectivity constraints to maintain network connectivity between agents.
     * 
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
template <typename T, unsigned int DIM>
class ConnectivityQPGenerator : public CBFQPGeneratorBase<T, DIM> {
  public:
    using typename CBFQPGeneratorBase<T, DIM>::VectorDIM;
    using typename CBFQPGeneratorBase<T, DIM>::Vector;
    using typename CBFQPGeneratorBase<T, DIM>::Row;
    using typename CBFQPGeneratorBase<T, DIM>::LinearConstraint;
    using typename CBFQPGeneratorBase<T, DIM>::Matrix;

    /**
         * @brief Constructor that takes a pointer to the ConnectivityCBF instance
         *
         * @param cbf Shared pointer to the Connectivity CBF implementation  
         * @param slack_config Configuration for different types of slack variables
         * @param num_neighbors Number of neighboring robots (for slack variable sizing)
         */
    explicit ConnectivityQPGenerator(std::shared_ptr<ConnectivityCBF> cbf,
                                     const SlackConfig& slack_config = SlackConfig{},
                                     size_t num_neighbors = 0);

    /**
         * @brief Adds connectivity constraint between agents to maintain network connectivity
         *
         * @param state Current state of the robot
         * @param robot_states Matrix of states for all robots in the network
         */
    void addConnConstraint(const Vector& state, const Eigen::MatrixXd& robot_states,
                           size_t self_idx);

    /**
         * @brief Adds CLF connectivity constraint
         *
         * @param state Current state of the robot
         * @param neighbor_state State of the target/neighboring robot
         * @param slack_idx Index of slack variable to use
         */
    void addCLFConstraint(const Vector& state, const Vector& neighbor_state, std::size_t slack_idx);

    /**
         * @brief Adds safety constraint between agents to prevent collisions
         *
         * @param state Current state of the robot
         * @param neighbor_state State of the target/neighboring robot
         * @param slack_idx Index of slack variable to use
         */
    void addSafetyConstraint(const Vector& state, const Vector& neighbor_state,
                             std::size_t slack_idx = 0);

    /**
         * @brief Adds minimum velocity constraints
         *
         * @param state Current state of the robot
         */
    void addMinVelConstraints(const Vector& state) override;

    /**
         * @brief Adds maximum velocity constraints
         *
         * @param state Current state of the robot
         */
    void addMaxVelConstraints(const Vector& state) override;

    // Slack configuration and variables
    SlackConfig slack_config_; ///< Configuration for different slack variable types
    std::vector<qpcpp::Variable<T>*>
        safety_slack_variables_; ///< Slack variables for safety constraints
    std::vector<qpcpp::Variable<T>*> clf_slack_variables_; ///< Slack variables for CLF constraints
    std::vector<qpcpp::Variable<T>*>
        connectivity_slack_variables_; ///< Slack variables for connectivity constraints

  private:
    std::shared_ptr<ConnectivityCBF> cbf_; ///< Pointer to the Connectivity CBF implementation

    // Using methods from the base class
    using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInput;
    using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInputWithSlackVariables;
};

} // namespace cbf

#endif //MPC_CONNECTIVITYQPGENERATOR_H