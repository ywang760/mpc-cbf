#ifndef CONNECTIVITY_CBF_H
#define CONNECTIVITY_CBF_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <math/Helpers.h>
#include <functional>
#include <ginac/ginac.h>
#include <vector>
#include <common/logging.hpp>

namespace cbf {

    class ConnectivityCBF {
    private:
        // Parameters
        double dmin, dmax;
        Eigen::VectorXd vmin, vmax;
        int STATE_VARS;
        int CONTROL_VARS;
        double gamma;
        double epsilon;
        // Symbols
        GiNaC::symbol px, py, th, vx, vy, w;                        // State for ego agent
        GiNaC::symbol px_n, py_n, vx_n, vy_n;                       // State for neighbor agent (omit th and w for simplicty)
        std::vector<GiNaC::symbol> px_list, py_list, eigenvec_list; // TODO: px_n and px_list contains duplicate information: room for optimization
        // System dynamics
        GiNaC::matrix A, B, f, g;
        // States for the ego and neighbor agent
        GiNaC::matrix state, p_n, v_n;
        // Symbolic constraints and bounds
        GiNaC::matrix Ac_safe, Ac_connectivity;
        GiNaC::matrix Ac_v1_max, Ac_v2_max, Ac_v3_max;
        GiNaC::matrix Ac_v1_min, Ac_v2_min, Ac_v3_min;
        GiNaC::ex Bc_safe, Bc_connectivity;
        GiNaC::ex Bc_v1_max, Bc_v2_max, Bc_v3_max;
        GiNaC::ex Bc_v1_min, Bc_v2_min, Bc_v3_min;
        // Alpha function
        std::function<GiNaC::ex(GiNaC::ex, double)> alpha;
        // Internal initialization
        std::pair<GiNaC::matrix, GiNaC::ex> initSafetyCBF();
        std::pair<GiNaC::matrix, GiNaC::ex> initVelCBF(GiNaC::ex bv);
        void initSymbolLists(int N);
        // Symbolic substitution utilities
        GiNaC::ex matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        GiNaC::ex valueSubs(GiNaC::ex m, Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        GiNaC::matrix matrixSubsMatrix(const GiNaC::matrix &expr_matrix,
                                       const Eigen::MatrixXd &robot_positions,
                                       const Eigen::VectorXd &eigenvec,
                                       const Eigen::Vector2d &self_position = Eigen::Vector2d::Zero());
        // Helpers for connectivity CBF
        GiNaC::matrix compute_dh_dx(int N, const GiNaC::ex& Rs, const GiNaC::ex& sigma);
        //GiNaC::matrix compute_d2h_dx2(const GiNaC::matrix& dh_dx_sym, int self_idx);
        Eigen::VectorXd compute_dLf_h_dx(
            const GiNaC::matrix& dh_dx_sym,
            int self_idx,
            const Eigen::MatrixXd& robot_positions,
            const Eigen::VectorXd& eigenvec,
            const Eigen::VectorXd& x_self,
            double Rs_val,
            double sigma_val);
        Eigen::Matrix2d compute_d2h_dx2_fd(
            const GiNaC::matrix& dh_dx_sym, 
            const Eigen::MatrixXd& robot_positions,
            const Eigen::VectorXd& eigenvec,
            const Eigen::Vector2d& x_self,
            int self_idx,
            double Rs_val,
            double sigma_val);

    public:
        using Vector3d = math::VectorDIM<double, 3>;
        ConnectivityCBF(double min_dist, double max_dist, Eigen::VectorXd vmin, Eigen::VectorXd vmax);
        ~ConnectivityCBF();
        // Basic constraints
        Eigen::VectorXd getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        double getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        double getMaxDistBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        // Velocity constraints
        Eigen::MatrixXd getMaxVelContraints(Eigen::VectorXd state);
        Eigen::MatrixXd getMinVelContraints(Eigen::VectorXd state);
        Eigen::VectorXd getMaxVelBounds(Eigen::VectorXd state);
        Eigen::VectorXd getMinVelBounds(Eigen::VectorXd state);
        std::pair<Eigen::VectorXd, double> initConnCBF(const Eigen::MatrixXd &robot_states,
                                                       const Eigen::VectorXd &x_self,
                                                       int self_idx);
        // Connectivity constraint //TODO: deprecated
        // Eigen::VectorXd getConnConstraints(const Eigen::VectorXd &x_self,
        //                                               const std::vector<Eigen::VectorXd> &other_positions);
        // double getConnBound(const Eigen::VectorXd &x_self,
        //                             const std::vector<Eigen::VectorXd> &other_positions);
        // Alpha setter
        void setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha);
    };
    // Free function
    std::pair<double, Eigen::VectorXd> getLambda2FromL(const Eigen::MatrixXd& robot_positions,
                                                       double Rs_value,
                                                       double sigma_value);
} // namespace cbf
#endif // CONNECTIVITY_CBF_H
