#ifndef CONNECTIVITY_CBF_H
#define CONNECTIVITY_CBF_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <functional>
#include <ginac/ginac.h>
#include <vector>
#include <common/logging.hpp>
#include <math/Types.h>   

namespace mpc_cbf {
    class ConnectivityCBF;
    GiNaC::ex     matrixSubs(GiNaC::matrix, Eigen::VectorXd, Eigen::VectorXd, const ConnectivityCBF&);
    GiNaC::ex     valueSubs (GiNaC::ex,     Eigen::VectorXd, Eigen::VectorXd, const ConnectivityCBF&);
    GiNaC::matrix matrixSubsMatrix(const GiNaC::matrix&,
                                       const Eigen::MatrixXd&,
                                       const Eigen::VectorXd&,
                                       const Eigen::Ref<const Eigen::VectorXd>&,
                                       const ConnectivityCBF&);
    } // namespace mpc_cbf

namespace mpc_cbf
{
template <typename T>
struct ConnectivityCBFParams {
    T dmin_, dmax_;
};
    class ConnectivityCBF
    {
        // 与 Helpers.hpp 中实现一致的 friend 声明（限定到 mpc_cbf::）
        friend GiNaC::ex     mpc_cbf::matrixSubs(GiNaC::matrix, Eigen::VectorXd, Eigen::VectorXd, const ConnectivityCBF&);
        friend GiNaC::ex     mpc_cbf::valueSubs (GiNaC::ex,     Eigen::VectorXd, Eigen::VectorXd, const ConnectivityCBF&);
        friend GiNaC::matrix mpc_cbf::matrixSubsMatrix(const GiNaC::matrix&,
                                                       const Eigen::MatrixXd&,
                                                       const Eigen::VectorXd&,
                                                       const Eigen::Ref<const Eigen::VectorXd>&,
                                                       const ConnectivityCBF&);


    private:
        // Parameters
        double dmin, dmax;
        Eigen::VectorXd vmin, vmax;
        int STATE_VARS;
        int CONTROL_VARS;
        double gamma;
        double epsilon;
        // Symbols
        GiNaC::symbol px, py, pz, vx, vy, vz;                        // State for ego agent
        GiNaC::symbol px_n, py_n, pz_n, vx_n, vy_n, vz_n;                       // State for neighbor agent (omit th and w for simplicty)
        std::vector<GiNaC::symbol> px_list, py_list, pz_list, eigenvec_list; // TODO: px_n and px_list contains duplicate information: room for optimization
        // System dynamics
        GiNaC::matrix A, B, f, g;
        // States for the ego and neighbor agent
        GiNaC::matrix state, p_n, v_n;
        // Symbolic constraints and bounds
        GiNaC::matrix Ac_safe, Ac_connectivity, Ac_CLF;
        GiNaC::matrix Ac_v1_max, Ac_v2_max, Ac_v3_max;
        GiNaC::matrix Ac_v1_min, Ac_v2_min, Ac_v3_min;
        GiNaC::ex Bc_safe, Bc_connectivity, Bc_CLF;
        GiNaC::ex Bc_v1_max, Bc_v2_max, Bc_v3_max;
        GiNaC::ex Bc_v1_min, Bc_v2_min, Bc_v3_min;
        // Alpha function
        std::function<GiNaC::ex(GiNaC::ex, double)> alpha;
        // Internal initialization
        std::pair<GiNaC::matrix, GiNaC::ex> initSafetyCBF();
        std::pair<GiNaC::matrix, GiNaC::ex> initCLFCBF();
        std::pair<GiNaC::matrix, GiNaC::ex> initVelCBF(GiNaC::ex bv);
        void initSymbolLists(int N);
        // Helpers for connectivity CBF
        GiNaC::matrix compute_dh_dx(int N, const GiNaC::ex &Rs, const GiNaC::ex &sigma);
        GiNaC::matrix compute_d2h_dx2(const GiNaC::matrix &dh_dx_sym, int self_idx);
        Eigen::VectorXd compute_dLf_h_dx(
            const GiNaC::matrix &dh_dx_sym,
            int self_idx,
            const Eigen::MatrixXd &robot_positions,
            const Eigen::VectorXd &eigenvec,
            const Eigen::VectorXd &x_self);
        // @quyichun check if the following functions are needed
        // Eigen::MatrixXd ginacToEigen(const GiNaC::matrix& m);
        // Eigen::Matrix2d compute_d2h_dx2_fd(
        //     const GiNaC::matrix& dh_dx_sym,
        //     const Eigen::MatrixXd& robot_positions,
        //     const Eigen::VectorXd& eigenvec,
        //     const Eigen::Vector2d& x_self,
        //     int self_idx,
        //     double Rs_val,
        //     double sigma_val);
        std::vector<std::vector<GiNaC::symbol>> Hij_masks_;

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
        double getEpsilon() const { return epsilon; }
        Eigen::VectorXd getCLFConstraints(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        double getCLFBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state);
        double getDmax() const { return dmax; }
        const std::vector<std::vector<GiNaC::symbol>>& HijMasks() const { return Hij_masks_; }
        void initMaskSymbols(int N);
    };
    // Free function
    std::pair<double, Eigen::VectorXd> getLambda2FromL(const Eigen::MatrixXd &robot_positions,
                                                       double Rs_value,
                                                       double sigma_value);
    
    
} // namespace mpc_cbf
#endif // CONNECTIVITY_CBF_H
