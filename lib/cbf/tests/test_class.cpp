#include <cbf/detail/cbf.h>


int main()
{
    double beta = 120 * M_PI / 180;
    double safety_dist = 3.0;
    double max_dist = 10.0;
    cbf::FovCBF cbf_test(beta, safety_dist, max_dist);

    Eigen::VectorXd state;
    state.resize(6);
    state.setZero();
    state(0) = -1;
    state(1) = -1;
    state(2) = 0.25*M_PI;


    Eigen::Vector2d target_state = 4*Eigen::Vector2d::Ones();
    target_state(1) = 2.5;

    Eigen::VectorXd Ac_safe = cbf_test.getSafetyConstraints(state, target_state);
    double Bc_safe = cbf_test.getSafetyBound(state, target_state);
    std::cout << "Ac_safe: " << Ac_safe.transpose() << std::endl;
    std::cout << "Bc_safe: " << Bc_safe << std::endl;

    Eigen::VectorXd Ac_lb = cbf_test.getLBConstraints(state, target_state);
    double Bc_lb = cbf_test.getLBBound(state, target_state);
    std::cout << "Ac_lb: " << Ac_lb.transpose() << std::endl;
    std::cout << "Bc_lb: " << Bc_lb << std::endl;

    Eigen::VectorXd Ac_rb = cbf_test.getRBConstraints(state, target_state);
    double Bc_rb = cbf_test.getRBBound(state, target_state);
    std::cout << "Ac_rb: " << Ac_rb.transpose() << std::endl;
    std::cout << "Bc_rb: " << Bc_rb << std::endl;

    Eigen::VectorXd Ac_range = cbf_test.getRangeConstraints(state, target_state);
    double Bc_range = cbf_test.getRangeBound(state, target_state);
    std::cout << "Ac_range: " << Ac_lb.transpose() << std::endl;
    std::cout << "Bc_range: " << Bc_lb << std::endl;




    return 0;

}