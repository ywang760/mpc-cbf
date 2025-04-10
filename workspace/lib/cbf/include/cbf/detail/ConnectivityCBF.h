#ifndef CONNECTIVITY_CBF_H
#define CONNECTIVITY_CBF_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <math/Helpers.h>
#include <functional>
#include <ginac/ginac.h>

namespace cbf
{
    template <typename T>
    struct ConnectivityCBFParams {
        T beta_;     // Convergence rate parameter
        T dmin_;     // Minimum connectivity distance
        T dmax_;     // Maximum connectivity distance
    };

    class ConnectivityCBF
    {
        private:
            // Constant params
            double dmin;         // Minimum connectivity distance
            double dmax;         // Maximum connectivity distance
            Eigen::VectorXd vmin;         // min velocity
            Eigen::VectorXd vmax;         // max velocity
            int STATE_VARS;
            int CONTROL_VARS;
            double gamma;

            GiNaC::symbol px, py, th, vx, vy, w, xt, yt;

            GiNaC::matrix state;
            GiNaC::matrix agent_state;
            GiNaC::matrix A;        // state matrix
            GiNaC::matrix B;        // input matrix
            GiNaC::matrix x;        // state variables
            GiNaC::matrix x_agent;  // agent variables
            GiNaC::matrix f;        // f(x)
            GiNaC::matrix g;        // g(x)
            std::vector<GiNaC::matrix> Ac_tot;              // Linear constraints
            std::vector<GiNaC::ex> Bc_tot;                  // upper bounds
            GiNaC::matrix Ac_safe;                          // Safety constraint
            GiNaC::matrix Ac_connectivity;                  // Connectivity constraint
            GiNaC::matrix Ac_v1_max;                        // Max velocity constraints
            GiNaC::matrix Ac_v2_max;
            GiNaC::matrix Ac_v3_max;
            GiNaC::matrix Ac_v1_min;                        // Min velocity constraints
            GiNaC::matrix Ac_v2_min;
            GiNaC::matrix Ac_v3_min;

            GiNaC::ex Bc_safe;                              // Safety bound
            GiNaC::ex Bc_connectivity;                      // Connectivity bound
            GiNaC::ex Bc_v1_max;                            // Max velocity bounds
            GiNaC::ex Bc_v2_max;
            GiNaC::ex Bc_v3_max;
            GiNaC::ex Bc_v1_min;                            // Min velocity bounds
            GiNaC::ex Bc_v2_min;
            GiNaC::ex Bc_v3_min;

            std::function<GiNaC::ex(GiNaC::ex, double)> alpha;

            std::pair<GiNaC::matrix, GiNaC::ex> initSafetyCBF();
            std::pair<GiNaC::matrix, GiNaC::ex> initConnectivityCBF();
            std::pair<GiNaC::matrix, GiNaC::ex> initVelCBF(GiNaC::ex bv);
            // TODO: these two are public helper functions -> could move to a helper class
            GiNaC::ex matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd agent_state);
            GiNaC::ex valueSubs(GiNaC::ex m, Eigen::VectorXd state, Eigen::VectorXd agent_state);

        public:
            ConnectivityCBF(double min_dist, double max_dist, Eigen::VectorXd vmin, Eigen::VectorXd vmax);
            ~ConnectivityCBF();
            Eigen::VectorXd getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd agent_state);
            Eigen::VectorXd getConnectivityConstraints(Eigen::VectorXd state, Eigen::VectorXd agent_state);
            Eigen::MatrixXd getMaxVelContraints(Eigen::VectorXd state);
            Eigen::MatrixXd getMinVelContraints(Eigen::VectorXd state);
            double getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd agent_state);
            double getMaxDistBound(Eigen::VectorXd state, Eigen::VectorXd agent_state);
            Eigen::VectorXd getMaxVelBounds(Eigen::VectorXd state);
            Eigen::VectorXd getMinVelBounds(Eigen::VectorXd state);
            void setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha);
    };
}

#endif // CONNECTIVITY_CBF_H