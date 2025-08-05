#include <iomanip>
#include <iostream>
#include <vector>

#include <particle_filter/detail/particle_filter.h>

namespace pf {
ParticleFilter::ParticleFilter() {
    std::cout << "Particle filter initialized.\n";
}

void ParticleFilter::init(int particles_num, Eigen::VectorXd initState, Eigen::MatrixXd initCov,
                          Eigen::MatrixXd processCov, Eigen::MatrixXd measurementCov) {
    n_ = particles_num;
    state_ = initState;
    W_ = processCov;     // process covariance for state update
    Q_ = initCov;        // actual particles distribution
    R_ = measurementCov; // measurement covariance
    dt = 0.2;
    std::srand(std::time(0)); // random seed for Eigen random gen

    // log_file_.open("samples.txt");

    state_size_ = state_.size();
    if (state_size_ != W_.cols() || state_size_ != W_.rows()) {
        throw std::runtime_error(
            "ParticleFilter Error: state size does not match state covariance size!");
    }
    if (state_size_ != Q_.cols() || state_size_ != Q_.rows()) {
        throw std::runtime_error(
            "ParticleFilter Error: state size does not match particles distribution size!");
    }

    // Create particles
    particles_.resize(state_size_, n_);

    // Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(Q_);
    Eigen::MatrixXd L = llt.matrixL();

    // Generate samples from std normal distrib
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < state_size_; ++j) {
            particles_(j, i) = dist(gen_);
        }
        particles_.col(i) = state_ + L * particles_.col(i);
    }

    // Init weights
    w_.resize(n_);
    w_.setOnes();
    w_ /= n_;

    std::cout << "Particle filter initialized\n";
}

ParticleFilter::~ParticleFilter() {
    std::cout << "Closing particle filter...\n";
}

// ---------------------- Particle Filter main loop -----------------------
void ParticleFilter::predict(Eigen::VectorXd input) {
    // Update particles according to state function
    particles_ = particles_.colwise() + input * dt;

    // Add process noise (Eigen random is uniform, not normal)
    // Eigen::MatrixXd noise = W_ * Eigen::MatrixXd::Random(state_size_, n_);
    // particles_ += noise;

    std::normal_distribution<double> dist(0.0, 1.0);
    auto func = [&](int) { return dist(gen_); };
    Eigen::MatrixXd noise = W_ * Eigen::MatrixXd::NullaryExpr(state_size_, n_, func);
    particles_ += noise;
}

void ParticleFilter::predict() {
    Eigen::VectorXd u;
    u.resize(state_size_);
    u.setZero();

    predict(u);
}

void ParticleFilter::update(Eigen::VectorXd measurement, Eigen::MatrixXd cov) {
    if (state_size_ != measurement.size()) {
        throw std::runtime_error(
            "ParticleFilter Error: state size does not match measurement size!");
    }
    if (state_size_ != cov.rows() || state_size_ != cov.cols()) {
        throw std::runtime_error(
            "ParticleFilter Error: state size does not match measurement covariance size!");
    }

    Eigen::MatrixXd innovation = particles_.colwise() - measurement;
    Eigen::MatrixXd dist = innovation.transpose() * cov.inverse() * innovation;
    Eigen::VectorXd coeff = dist.diagonal();
    w_ = (-0.5 * coeff).array().exp();
    w_ /= w_.sum();
}

void ParticleFilter::update(Eigen::VectorXd measurement) {
    update(measurement, R_);
}

void ParticleFilter::resample() {
    std::vector<double> weights(w_.data(), w_.data() + w_.size());
    Eigen::MatrixXd resampled_particles;
    resampled_particles.resize(state_size_, n_);
    std::discrete_distribution<int> distribution(weights.begin(), weights.end());
    for (int i = 0; i < n_; ++i) {
        int index = distribution(gen_);
        resampled_particles.col(i) = particles_.col(index);
        w_(i) = weights[index];
    }

    particles_ = resampled_particles;
}

void ParticleFilter::estimateState() {
    // Update state and distribution
    state_ = particles_.rowwise().mean();
    Q_ = getDistribution();
}

// ----------------- End of Particle Filter main loop ------------

void ParticleFilter::saveParticles() {
    if (log_file_.is_open()) {
        for (int i = 0; i < particles_.cols(); ++i) {
            for (int j = 0; j < particles_.rows(); ++j) {
                log_file_ << std::setw(5) << std::setfill('0') << particles_(j, i) << ',';
            }
            log_file_ << "\n";
        }
    }
}

void ParticleFilter::setParticles(Eigen::MatrixXd particles) {
    if (particles.rows() != state_size_) {
        throw std::runtime_error("ParticleFilter Error: Dimension mismatch!");
    }
    if (particles.cols() != n_) {
        throw std::runtime_error("ParticleFilter Error: Dimension mismatch!");
    }
    particles_ = particles;
}

Eigen::MatrixXd ParticleFilter::getParticles() {
    return particles_;
}

Eigen::MatrixXd ParticleFilter::getDistribution() {
    Eigen::MatrixXd cov_matrix(state_size_, state_size_);
    cov_matrix.setZero();

    /*
        for (int i = 0; i < n_; ++i)
        {
            Eigen::Vector3d p = particles_.col(i);
            Eigen::Vector3d diff = p - state_;
            cov_matrix += diff * diff.transpose();
        }
        */

    Eigen::MatrixXd diff = particles_.colwise() - state_;
    cov_matrix += diff * diff.transpose();
    cov_matrix /= (n_ - 1);

    return cov_matrix;
}

Eigen::VectorXd ParticleFilter::getState() {
    return state_;
}

Eigen::VectorXd ParticleFilter::getWeights() {
    return w_;
}

void ParticleFilter::setWeights(Eigen::VectorXd weights) {
    if (weights.size() != n_) {
        throw std::runtime_error(
            "ParticleFilter Error: Weights size different from number of particles!");
    }
    w_ = weights;
}

void ParticleFilter::setWeights(std::vector<double> weights) {
    if (weights.size() != n_) {
        throw std::runtime_error(
            "ParticleFilter Error: Weights size different from number of particles!");
    }
    w_ = Eigen::VectorXd::Map(&weights[0], weights.size());
}

} // namespace pf