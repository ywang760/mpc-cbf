#include <math/Random.h>
#include <random>

namespace math {
    // Template implementation
    template <typename T, unsigned int DIM>
    model::State<T, DIM> addRandomNoise(const model::State<T, DIM> &state, T pos_std, T vel_std)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<T> distribution_position(T(0.0), pos_std);
        std::normal_distribution<T> distribution_velocity(T(0.0), vel_std);

        model::State<T, DIM> result_state = state;

        // Add noise to position components
        for (unsigned int i = 0; i < DIM; ++i)
        {
            T pos_sample = distribution_position(gen);
            result_state.pos_(i) += pos_sample;
        }

        // Add noise to velocity components
        for (unsigned int i = 0; i < DIM; ++i)
        {
            T vel_sample = distribution_velocity(gen);
            result_state.vel_(i) += vel_sample;
        }

        return result_state;
    }

    // Explicit instantiations
    template model::State<double, 2U> addRandomNoise<double, 2U>(
        const model::State<double, 2U> &state, double pos_std, double vel_std);
    template model::State<double, 3U> addRandomNoise<double, 3U>(
        const model::State<double, 3U> &state, double pos_std, double vel_std);
    template model::State<float, 2U> addRandomNoise<float, 2U>(
        const model::State<float, 2U> &state, float pos_std, float vel_std);
    template model::State<float, 3U> addRandomNoise<float, 3U>(
        const model::State<float, 3U> &state, float pos_std, float vel_std);

} // namespace math