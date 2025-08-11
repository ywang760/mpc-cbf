//
// Created on 4/10/25.
//

#ifndef CBF_CONNECTIVITYCONTROL_H
#define CBF_CONNECTIVITYCONTROL_H

#include <cbf/detail/ConnectivityCBF.h>
#include <cbf/optimization/ConnectivityQPGenerator.h>
#include <model/DoubleIntegrator.h>
#include <numeric>
#include <qpcpp/solvers/CPLEX.h>

namespace cbf {
template <typename T, unsigned int DIM>
class ConnectivityControl {
  public:
    using QPGenerator = cbf::ConnectivityQPGenerator<T, DIM>;
    using Problem = qpcpp::Problem<T>;
    using CPLEXSolver = qpcpp::CPLEXSolver<T>;
    using SolveStatus = qpcpp::SolveStatus;
    using Vector = math::Vector<T>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using Matrix = math::Matrix<T>;
    using State = model::State<T, DIM>;

    /**
         * @brief Constructor for ConnectivityControl
         * 
         * @param cbf Shared pointer to ConnectivityCBF implementation
         * @param slack_config Configuration for different types of slack variables
         * @param num_neighbors Number of neighboring robots (for slack variable sizing)
         */
    ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf, 
                        const SlackConfig& slack_config = SlackConfig{}, 
                        size_t num_neighbors = 0);
    ~ConnectivityControl() = default;

    /**
         * @brief Solve the optimization problem to get connectivity-preserving control input
         *
         * @param cbf_u Output parameter for the optimized control input
         * @param desired_u Desired control input without CBF constraints
         * @param current_states Current states of all robots in the network
         * @param self_idx Index of the ego robot for which we are computing control input
         * @param u_min Minimum allowable control inputs
         * @param u_max Maximum allowable control inputs
         * @return true if optimization was successful
         * @return false if optimization failed
         */
    bool optimize(VectorDIM& cbf_u, const VectorDIM& desired_u, std::vector<State> current_states,
                  size_t self_idx, const VectorDIM& u_min, const VectorDIM& u_max);

  private:
    QPGenerator qp_generator_;
    SlackConfig slack_config_;
    size_t num_neighbors_;
    std::shared_ptr<ConnectivityCBF> cbf_;
};

} // namespace cbf

#endif //CBF_CONNECTIVITYCONTROL_H