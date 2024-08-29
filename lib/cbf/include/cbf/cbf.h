
#ifndef FOV_CBF_H
#define FOV_CBF_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <ginac/ginac.h>


namespace cbf
{   
    // template <typename T>
    class FovCBF
    {
        private:
            // Constant params
            double fov;          // FoV
            double Ds;           // Safety dist
            double Rs;           // Max dist
            int STATE_VARS;
            int CONTROL_VARS;

            GiNaC::symbol px, py, th, vx, vy, w, xt, yt;
            
            GiNaC::matrix state;
            GiNaC::matrix target_state;
            GiNaC::matrix A;        // state matrix
            GiNaC::matrix B;        // input matrix
            GiNaC::matrix x;        // state variables
            GiNaC::matrix x_target; // target variables
            GiNaC::matrix f;        // f(x)
            GiNaC::matrix g;        // g(x)
            std::vector<GiNaC::matrix> Ac_tot;              // Linear constraints
            std::vector<GiNaC::ex> Bc_tot;              // upper bounds
            GiNaC::matrix Ac_safe;
            GiNaC::matrix Ac_lb;                        // left border constraints
            GiNaC::matrix Ac_rb;                        // right border constraints
            GiNaC::matrix Ac_range;
            GiNaC::ex Bc_safe;
            GiNaC::ex Bc_lb;                            // left border
            GiNaC::ex Bc_rb;                            // right border
            GiNaC::ex Bc_range;

            std::pair<GiNaC::matrix, GiNaC::ex> initSafetyCBF();
            std::pair<GiNaC::matrix, GiNaC::ex> initBorder1CBF();
            std::pair<GiNaC::matrix, GiNaC::ex> initBorder2CBF();
            std::pair<GiNaC::matrix, GiNaC::ex> initRangeCBF();
            GiNaC::ex matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd target_state);
            GiNaC::ex valueSubs(GiNaC::ex m, Eigen::VectorXd state, Eigen::VectorXd target_state);


        public:
            FovCBF(double fov, double safety_dist, double max_dist);
            ~FovCBF();
            Eigen::VectorXd getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state);
            Eigen::VectorXd getRangeConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state);
            Eigen::VectorXd getLBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state);
            Eigen::VectorXd getRBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state);
            double getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd target_state);
            double getRangeBound(Eigen::VectorXd state, Eigen::VectorXd target_state);
            double getLBBound(Eigen::VectorXd state, Eigen::VectorXd target_state);
            double getRBBound(Eigen::VectorXd state, Eigen::VectorXd target_state);



    };
}

#endif // FOV_CBF_H