#ifndef MPC_FOVQPGENERATOR_H
#define MPC_FOVQPGENERATOR_H

#include <cbf/detail/FovCBF.h>
#include <cbf/optimization/CBFQPGeneratorBase.h>

namespace cbf {
/**
     * @brief Field of View Quadratic Program Generator
     * 
     * This class implements the FOV-specific constraint methods for CBF QP generation.
     * It handles field-of-view constraints like left/right borders and range constraints.
     * 
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
template <typename T, unsigned int DIM>
class FovQPGenerator : public CBFQPGeneratorBase<T, DIM> {
  public:
    using typename CBFQPGeneratorBase<T, DIM>::VectorDIM;
    using typename CBFQPGeneratorBase<T, DIM>::Vector;
    using typename CBFQPGeneratorBase<T, DIM>::Row;
    using typename CBFQPGeneratorBase<T, DIM>::LinearConstraint;
    using typename CBFQPGeneratorBase<T, DIM>::Matrix;

    /**
         * @brief Constructor that takes a pointer to the FovCBF instance
         *
         * @param cbf Shared pointer to the Field of View CBF implementation
         * @param num_robots Number of neighboring agents (used for slack variables)
         * @param slack_mode If true, add slack variables to allow constraint relaxation
         */
    explicit FovQPGenerator(std::shared_ptr<FovCBF> cbf, int num_robots = 0,
                            bool slack_mode = false);

    /**
         * @brief Adds safety constraint between agents to prevent collisions
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
    void addSafetyConstraint(const Vector& state, const Vector& target_state,
                             bool use_slack = false, std::size_t slack_idx = 0) override;

    /**
         * @brief Adds field-of-view left border constraint
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
    void addLeftBorderConstraint(const Vector& state, const Vector& target_state,
                                 bool use_slack = false, std::size_t slack_idx = 0);

    /**
         * @brief Adds field-of-view right border constraint
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
    void addRightBorderConstraint(const Vector& state, const Vector& target_state,
                                  bool use_slack = false, std::size_t slack_idx = 0);

    /**
         * @brief Adds range constraint to maintain target within sensor range
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
    void addRangeConstraint(const Vector& state, const Vector& target_state, bool use_slack = false,
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

  private:
    std::shared_ptr<FovCBF> cbf_; ///< Pointer to the Field of View CBF implementation

    // Using methods from the base class
    using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInput;
    using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInputWithSlackVariables;
};

} // namespace cbf

#endif //MPC_FOVQPGENERATOR_H