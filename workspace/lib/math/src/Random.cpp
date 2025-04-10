#include <math/Random.h>
#include <random>

namespace math {

Vector<double> addRandomNoise(const Vector<double>& xt, double pos_std, double vel_std, unsigned int dim) {
    std::random_device rd;
    std::mt19937 gen(rd());
    Vector<double> result_xt = xt;
    double sample = 0.0;
    std::normal_distribution<double> distribution_position(0.0, pos_std);
    std::normal_distribution<double> distribution_velocity(0.0, vel_std);
    for (int i = 0; i < 2*dim; ++i) {
        if (i < dim) {
            sample = distribution_position(gen);
        } else {
            sample = distribution_velocity(gen);
        }
        result_xt(i) += sample;
    }
    return result_xt;
}

} // namespace math
