#include <cbf/detail/cbf.h>


int main()
{
    double beta = 120 * M_PI / 180;
    double safety_dist = 2.0;
    double max_dist = 10.0;
    cbf::FovCBF cbf_test(beta, safety_dist, max_dist);
    
    Eigen::VectorXd state;
    state.resize(6);
    state.setZero();

    Eigen::Vector2d target_state = 4*Eigen::Vector2d::Ones();

    Eigen::VectorXd Ac = cbf_test.getLBConstraints(state, target_state);
    double Bc = cbf_test.getLBBound(state, target_state);

    std::cout << "Ac: " << Ac.transpose() << std::endl;
    std::cout << "Bc: " << Bc << std::endl;


    return 0;

}