#include <cbf/detail/FovCBF.h>

int main()
{
    double beta = 120 * M_PI / 180;
    double safety_dist = 3.0;
    double max_dist = 10.0;
    double vmax = 2.0;
    cbf::FovCBF cbf_test(beta, safety_dist, max_dist, vmax);
    
    Eigen::VectorXd state;
    state.resize(6);
    state.setZero();
    state(0) = -1;
    state(1) = -1;
    state(2) = 0.25 * M_PI;
    state(3) = 1.0;

    Eigen::Vector2d target_state = 4*Eigen::Vector2d::Ones();
    target_state(1) = 2.5;

    Eigen::VectorXd Ac_safe = cbf_test.getSafetyConstraints(state, target_state);
    double Bc_safe = cbf_test.getSafetyBound(state, target_state);

    Eigen::VectorXd Ac_lb = cbf_test.getLBConstraints(state, target_state);
    double Bc_lb = cbf_test.getLBBound(state, target_state);

    Eigen::VectorXd Ac_rb = cbf_test.getRBConstraints(state, target_state);
    double Bc_rb = cbf_test.getRBBound(state, target_state);

    Eigen::VectorXd Ac_range = cbf_test.getRangeConstraints(state, target_state);
    double Bc_range = cbf_test.getRangeBound(state, target_state);

    Eigen::MatrixXd Ac_vmax = cbf_test.getMaxVelContraints(state);
    Eigen::MatrixXd Ac_vmin = cbf_test.getMinVelContraints(state);
    Eigen::VectorXd Bc_vmax = cbf_test.getMaxVelBounds(state);
    Eigen::VectorXd Bc_vmin = cbf_test.getMinVelBounds(state);

    std::cout << "Max vel matrix: \n" << Ac_vmax << std::endl;
    std::cout << "Min vel matrix: \n" << Ac_vmin << std::endl;
    std::cout << "Max vel bounds: " << Bc_vmax.transpose() << std::endl;
    std::cout << "Min vel bounds: " << Bc_vmin.transpose() << std::endl;




    return 0;

}