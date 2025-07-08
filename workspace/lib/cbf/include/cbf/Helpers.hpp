#ifndef CBF_HELPERS_H
#define CBF_HELPERS_H

#include <eigen3/Eigen/Dense>
#include <ginac/ginac.h>
#include <cbf/detail/ConnectivityCBF.h>

namespace cbf {

    inline GiNaC::ex matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd neighbor_state,
                                 const ConnectivityCBF& cbf) {
        // Substitute each state variable with its numerical value
        GiNaC::ex tmp = GiNaC::subs(a, cbf.px == state(0));
        tmp = GiNaC::subs(tmp, cbf.py == state(1));
        tmp = GiNaC::subs(tmp, cbf.th == state(2));
        tmp = GiNaC::subs(tmp, cbf.vx == state(3));
        tmp = GiNaC::subs(tmp, cbf.vy == state(4));
        tmp = GiNaC::subs(tmp, cbf.w == state(5));
        tmp = GiNaC::subs(tmp, cbf.px_n == neighbor_state(0));
        tmp = GiNaC::subs(tmp, cbf.py_n == neighbor_state(1));
        // pz_n (th_n) is omitted for now
        tmp = GiNaC::subs(tmp, cbf.vx_n == neighbor_state(3));
        tmp = GiNaC::subs(tmp, cbf.vy_n == neighbor_state(4));
        // w_n (angular velocity of neighbor) is omitted for now
        return tmp;
    }

    inline GiNaC::ex valueSubs(GiNaC::ex a, Eigen::VectorXd state, Eigen::VectorXd neighbor_state,
                               const ConnectivityCBF& cbf) {
        // Substitute each state variable with its numerical value
        GiNaC::ex tmp = GiNaC::subs(a, cbf.px == state(0));
        tmp = GiNaC::subs(tmp, cbf.py == state(1));
        tmp = GiNaC::subs(tmp, cbf.th == state(2));
        tmp = GiNaC::subs(tmp, cbf.vx == state(3));
        tmp = GiNaC::subs(tmp, cbf.vy == state(4));
        tmp = GiNaC::subs(tmp, cbf.w == state(5));
        tmp = GiNaC::subs(tmp, cbf.px_n == neighbor_state(0));
        tmp = GiNaC::subs(tmp, cbf.py_n == neighbor_state(1));
        // pz_n (th_n) is omitted for now
        tmp = GiNaC::subs(tmp, cbf.vx_n == neighbor_state(3));
        tmp = GiNaC::subs(tmp, cbf.vy_n == neighbor_state(4));
        // w_n (angular velocity of neighbor) is omitted for now
        return tmp;
    }

    inline GiNaC::matrix matrixSubsMatrix(
        const GiNaC::matrix &expr_matrix,
        const Eigen::MatrixXd &robot_positions,
        const Eigen::VectorXd &eigenvec,
        const Eigen::Vector2d &self_position,
        const ConnectivityCBF &cbf)
    {
        GiNaC::exmap substitutions;
        // 机器人位置和特征值向量替换
        for (int i = 0; i < robot_positions.rows(); ++i)
        {
            substitutions[cbf.px_list[i]] = robot_positions(i, 0);
            substitutions[cbf.py_list[i]] = robot_positions(i, 1);
            substitutions[cbf.eigenvec_list[i]] = eigenvec(i);
        }
        // 当前机器人的自身位置
        substitutions[cbf.px] = self_position(0);
        substitutions[cbf.py] = self_position(1);
        // 对整个矩阵统一替换
        GiNaC::matrix result = GiNaC::ex_to<GiNaC::matrix>(expr_matrix.subs(substitutions));
        return result;
    }
} // namespace cbf


#endif // CBF_HELPERS_H