#include <iostream>
#include <fstream>


#include <particle_filter/detail/particle_filter.h>

int main()
{
    int samples_num = 100;
    Eigen::VectorXd state = Eigen::Vector3d::Ones();
    std::cout << "Init state: " << state << std::endl;
    Eigen::MatrixXd initCov = 1.0*Eigen::Matrix3d::Identity();
    Eigen::MatrixXd processCov = 1.0*Eigen::Matrix3d::Identity();
    Eigen::MatrixXd measCov = 0.5*Eigen::Matrix3d::Identity();

    pf::ParticleFilter filter;
    filter.init(samples_num, state, initCov, processCov, measCov);
    std::cout << "Dsitrib: \n" << filter.getDistribution() << std::endl;
    // Eigen::Vector3d detection;
    // detection << 2.3, 1.5, 1.1;
    std::random_device r;
    std::default_random_engine gen;

    for (int i = 0; i < 10; ++i)
    { 
        double sigma = 0.5;
        std::normal_distribution<double> dx(state(0), sigma);
        std::normal_distribution<double> dy(state(1), sigma);
        std::normal_distribution<double> dz(state(2), sigma);
        Eigen::Vector3d detection;
        detection << dx(gen), dy(gen), dz(gen);   
        std::cout << "detection at t = " << i << ": " << detection.transpose() << std::endl;
        filter.predict();
        // filter.update(detection);
        filter.resample();
        filter.estimateState();
        filter.saveParticles();
        std::cout << "Estimate at t = " << i << ": " << filter.getState().transpose() << std::endl;

    }
}