//
// Created by lishuo on 1/31/24.
//

#include <splines/curves/Bezier.h>

namespace splines {
    template<typename T, unsigned int DIM>
    Bezier<T, DIM>::Bezier(T max_parameter, const std::vector<VectorDIM> &control_points)
            : max_parameter_(max_parameter), control_points_(control_points) {}

    template<typename T, unsigned int DIM>
    void Bezier<T, DIM>::appendControlPoint(const VectorDIM &control_point) {
        control_points_.push_back(control_point);
    }

    template<typename T, unsigned int DIM>
    std::size_t Bezier<T, DIM>::numControlPoints() const {
        return control_points_.size();
    }

    template<typename T, unsigned int DIM>
    const typename Bezier<T, DIM>::VectorDIM
    Bezier<T, DIM>::getControlPoint(std::size_t control_point_idx) const {
        return control_points_[control_point_idx];
    }

    template<typename T, unsigned int DIM>
    const typename std::vector<typename Bezier<T, DIM>::VectorDIM>
    Bezier<T, DIM>::getControlPoints() const {
        std::vector<Bezier<T, DIM>::VectorDIM> result;

        for (const VectorDIM &control_point: control_points_) {
            result.push_back(control_point);
        }

        return result;
    }

    template<typename T, unsigned int DIM>
    uint64_t Bezier<T, DIM>::degree() const {
        return control_points_.size() - 1;
    }

    template<typename T, unsigned int DIM>
    T Bezier<T, DIM>::max_parameter() const {
        return max_parameter_;
    }

    template<typename T, unsigned int DIM>
    void Bezier<T, DIM>::set_max_parameter(T max_parameter) {
        max_parameter_ = max_parameter;
    }

    template<typename T, unsigned int DIM>
    void Bezier<T, DIM>::scaleMaxParameter(T scaling_factor) {
        max_parameter_ *= scaling_factor;
    }

    template<typename T, unsigned int DIM>
    typename Bezier<T, DIM>::VectorDIM
    Bezier<T, DIM>::eval(T parameter, uint64_t derivative_degree) const {
        using Row = math::Row<T>;
        // construct the bernsteinBasis
        Row basis_functions = bernsteinBasis(control_points_.size() - 1, max_parameter_,
                                              parameter, derivative_degree);

        VectorDIM result = VectorDIM::Zero();
        for (std::size_t control_point_idx = 0; control_point_idx < control_points_.size(); ++control_point_idx) {
            result += control_points_[control_point_idx] * basis_functions(control_point_idx);
        }

        return result;
    }

    template<typename T, unsigned int DIM>
    typename Bezier<T, DIM>::MaximumDerivativeMagnitude
    Bezier<T, DIM>::maximumDerivativeMagnitude(uint64_t derivative_degree) const {
        constexpr T step_size = 0.01;
        T maximum_magnitude = std::numeric_limits<T>::lowest();
        T maximum_parameter = 0;
        for (T parameter = 0; parameter <= max_parameter_; parameter += step_size) {
            const VectorDIM &eval_result = eval(parameter, derivative_degree);

            const T eval_result_norm = eval_result.norm();
            if (eval_result_norm > maximum_magnitude) {
                maximum_parameter = parameter;
                maximum_magnitude = eval_result_norm;
            }
        }

        return MaximumDerivativeMagnitude{
                .magnitude = maximum_magnitude,
                .parameter = maximum_parameter,
        };
    }

    template
    class Bezier<double, 2U>;

    template
    class Bezier<float, 2U>;

    template
    class Bezier<double, 3U>;

    template
    class Bezier<float, 3U>;

} // splines