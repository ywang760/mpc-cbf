#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <random>

namespace pf {
class ParticleFilter {
  private:
    int n_;                     // number of
    size_t state_size_;         // state size
    Eigen::VectorXd state_;     // state of the robot
    Eigen::VectorXd w_;         // samples weights
    Eigen::MatrixXd particles_; // samples matrix: 1 particle (x,y,th) per column
    Eigen::MatrixXd W_;         // process noise
    Eigen::MatrixXd Q_;         // samples distribution
    Eigen::MatrixXd R_;         // measurement noise

    double dt;

    std::mt19937 gen_; // random generator
    std::ofstream log_file_;

  public:
    ParticleFilter();
    ~ParticleFilter();
    void init(int particles_num, Eigen::VectorXd initState, Eigen::MatrixXd initCov,
              Eigen::MatrixXd processCov, Eigen::MatrixXd measurementCov);
    void predict();
    void predict(Eigen::VectorXd input);
    void update(Eigen::VectorXd observation, Eigen::MatrixXd measurementCov);
    void update(Eigen::VectorXd observation);
    void resample();
    void estimateState();
    Eigen::MatrixXd getParticles();
    void setParticles(Eigen::MatrixXd particles);
    Eigen::VectorXd getWeights();
    Eigen::VectorXd getState();
    Eigen::MatrixXd getDistribution();
    void saveParticles();
    void setWeights(Eigen::VectorXd weights);
    void setWeights(std::vector<double> weights);
};
} // namespace pf

#endif //PARTICLE_FILTER_H