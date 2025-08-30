#ifndef CONNECTIVITY_CBF_H
#define CONNECTIVITY_CBF_H

#include <common/logging.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <functional>
#include <ginac/ginac.h>
#include <math/Helpers.h>
#include <vector>

namespace cbf {

/**
 * @brief Configuration structure for slack variables in connectivity CBF constraints
 */
struct SlackConfig {
    bool safety_slack = false;              // Enable slack for safety constraints
    bool clf_slack = false;                 // Enable slack for CLF constraints  
    bool connectivity_slack = false;        // Enable slack for connectivity constraints
    double safety_slack_cost = 100000.0;    // Cost penalty for safety constraint violations
    double clf_slack_cost = 50000.0;        // Cost penalty for CLF constraint violations
    double connectivity_slack_cost = 25000.0; // Cost penalty for connectivity constraint violations
    double slack_decay_rate = 0.1;            // Slack decay rate over prediction horizon
};

template <typename T>
struct ConnectivityCBFParams {
    T dmin_, dmax_;
    SlackConfig slack_config_;
};

class ConnectivityCBF {
    // Friend declarations for helper functions
    friend GiNaC::ex matrixSubs(GiNaC::matrix a, Eigen::VectorXd state,
                                Eigen::VectorXd neighbor_state, const ConnectivityCBF& cbf);
    friend GiNaC::ex valueSubs(GiNaC::ex a, Eigen::VectorXd state, Eigen::VectorXd neighbor_state,
                               const ConnectivityCBF& cbf);
    friend GiNaC::matrix matrixSubsFull(const GiNaC::matrix& expr_matrix,
                                        const Eigen::MatrixXd& robot_states,
                                        const Eigen::VectorXd& eigenvec,
                                        const Eigen::VectorXd& state, const ConnectivityCBF& cbf);
    friend GiNaC::ex valueSubsFull(const GiNaC::ex& expr, const Eigen::MatrixXd& robot_states,
                                   const Eigen::VectorXd& eigenvec, const Eigen::VectorXd& state,
                                   const ConnectivityCBF& cbf);

  private:
    // Parameters
    double dmin, dmax;
    Eigen::VectorXd vmin, vmax;
    int STATE_VARS;
    int CONTROL_VARS;
    double gamma;
    double epsilon;
    // Symbols
    GiNaC::symbol px, py, th, vx, vy, w;  // State for ego agent
    GiNaC::symbol px_n, py_n, vx_n, vy_n; // State for neighbor agent (omit th and w for simplicty)
    std::vector<GiNaC::symbol> px_list, py_list,
        eigenvec_list; // TODO: px_n and px_list contains duplicate information: room for optimization
    // System dynamics
    GiNaC::matrix A, B, f, g;
    // States for the ego and neighbor agent
    GiNaC::matrix state, p_n, v_n;
    // Symbolic constraints and bounds
    GiNaC::matrix Ac_safe, Ac_conn, Ac_CLF;
    GiNaC::matrix Ac_v1_max, Ac_v2_max, Ac_v3_max;
    GiNaC::matrix Ac_v1_min, Ac_v2_min, Ac_v3_min;
    GiNaC::ex Bc_safe, Bc_conn, Bc_CLF;
    GiNaC::ex Bc_v1_max, Bc_v2_max, Bc_v3_max;
    GiNaC::ex Bc_v1_min, Bc_v2_min, Bc_v3_min;

    GiNaC::symbol h;

    // Alpha function
    std::function<GiNaC::ex(GiNaC::ex, double)> alpha;
    // Internal initialization
    std::pair<GiNaC::matrix, GiNaC::ex> initSafetyCBF();
    std::pair<GiNaC::matrix, GiNaC::ex> initCLFCBF();
    std::pair<GiNaC::matrix, GiNaC::ex> initVelCBF(GiNaC::ex bv);
    void initSymbolLists(int N);
    // Helpers for connectivity CBF
    GiNaC::matrix compute_full_grad_h(int N, const GiNaC::ex& Rs, const GiNaC::ex& sigma);

  public:
    using Vector3d = math::VectorDIM<double, 3>;
    ConnectivityCBF(double min_dist, double max_dist, Eigen::VectorXd vmin, Eigen::VectorXd vmax);
    ~ConnectivityCBF();
    // Safety constraints
    Eigen::VectorXd getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
    double getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
    // Connectivity Constraints (CBF)
    double getSigma() const;
    std::pair<double, Eigen::VectorXd> getLambda2(const Eigen::MatrixXd& robot_positions);
    Eigen::VectorXd getConnConstraints(Eigen::VectorXd state, Eigen::MatrixXd robot_states,
                                       Eigen::VectorXd eigenvec);
    double getConnBound(Eigen::VectorXd state, Eigen::MatrixXd robot_states,
                        Eigen::VectorXd eigenvec, double h_val);
    // Connectivity Constraints (CLF)
    Eigen::VectorXd getCLFConstraints(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
    double getCLFBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
    // Velocity constraints
    Eigen::MatrixXd getMaxVelContraints(Eigen::VectorXd state);
    Eigen::MatrixXd getMinVelContraints(Eigen::VectorXd state);
    Eigen::VectorXd getMaxVelBounds(Eigen::VectorXd state);
    Eigen::VectorXd getMinVelBounds(Eigen::VectorXd state);
    // TODO: initConnCBF should be private
    std::pair<GiNaC::matrix, GiNaC::ex> initConnCBF(int N, int self_idx);
    // Alpha setter
    void setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha);
};
// Free function
} // namespace cbf
#endif // CONNECTIVITY_CBF_H
