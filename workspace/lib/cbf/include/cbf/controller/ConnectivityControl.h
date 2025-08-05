//
// Created on 4/10/25.
//

#ifndef CBF_CONNECTIVITYCONTROL_H
#define CBF_CONNECTIVITYCONTROL_H

#include <cbf/detail/ConnectivityCBF.h>
#include <cbf/optimization/ConnectivityQPGenerator.h>
#include <model/DoubleIntegrator.h>
#include <qpcpp/solvers/CPLEX.h>
#include <numeric>

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
         * @param number_neighbors Number of neighboring agents to consider
         * @param slack_mode Whether to use slack variables in constraints
         * @param slack_cost Cost coefficient for slack variables
         * @param slack_decay_rate Decay rate for slack variables
         */
        ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf, int number_neighbors=0, bool slack_mode=false, T slack_cost=1000, T slack_decay_rate=0.2);
        ~ConnectivityControl()=default;

        /**
         * @brief Solve the optimization problem to get connectivity-preserving control input
         *
         * @param cbf_u Output parameter for the optimized control input
         * @param desired_u Desired control input without CBF constraints
         * @param current_states Current states of all robots in the network
         * @param ego_robot_idx Index of the ego robot for which we are computing control input
         * @param u_min Minimum allowable control inputs
         * @param u_max Maximum allowable control inputs
         * @return true if optimization was successful
         * @return false if optimization failed
         */
        bool optimize(VectorDIM &cbf_u,
                      const VectorDIM &desired_u,
                      std::vector<State> current_states,
                      size_t ego_robot_idx,
                      const VectorDIM &u_min,
                      const VectorDIM &u_max);

        const std::shared_ptr<cbf::ConnectivityCBF>& getCBF() const {
            return qp_generator_.getCBF();
        }


    private:
        QPGenerator qp_generator_;
        bool slack_mode_;
        T slack_cost_;
        T slack_decay_rate_;
        std::shared_ptr<ConnectivityCBF> cbf_;
    };

} // cbf

#endif //CBF_CONNECTIVITYCONTROL_H