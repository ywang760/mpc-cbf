//
// Created by lishuo on 2/13/24.
//

#include <splines/curves/SingleParameterPiecewiseCurve.h>

namespace splines {

template <typename T, unsigned int DIM>
void SingleParameterPiecewiseCurve<T, DIM>::addPiece(
    std::unique_ptr<SingleParameterCurve>&& piece_ptr) {
    T max_parameter_of_piece = piece_ptr->max_parameter();
    if (cumulative_max_parameters_.empty()) {
        cumulative_max_parameters_.push_back(max_parameter_of_piece);
    } else {
        cumulative_max_parameters_.push_back(cumulative_max_parameters_.back() +
                                             max_parameter_of_piece);
    }
    piece_ptrs_.emplace_back(std::move(piece_ptr));
}

template <typename T, unsigned int DIM>
std::size_t SingleParameterPiecewiseCurve<T, DIM>::numPieces() const {
    return piece_ptrs_.size();
}

template <typename T, unsigned int DIM>
std::reference_wrapper<const SingleParameterCurve<T, DIM>>
SingleParameterPiecewiseCurve<T, DIM>::getPiece(std::size_t piece_idx) const {
    if (piece_idx >= piece_ptrs_.size()) {
        throw std::out_of_range(
            "SingleParameterPiecewiseCurve::getPiece: piece index out of range.");
    }

    SingleParameterCurve* piece_ptr = piece_ptrs_[piece_idx].get();
    return *piece_ptr;
}

template <typename T, unsigned int DIM>
void SingleParameterPiecewiseCurve<T, DIM>::setMaxParameter(std::size_t piece_idx,
                                                            T new_max_parameter) {
    if (piece_idx >= piece_ptrs_.size()) {
        throw std::out_of_range(
            "SingleParameterPiecewiseCurve::setMaxParameter: piece index out of range.");
    }

    piece_ptrs_[piece_idx]->set_max_parameter(new_max_parameter);

    for (std::size_t idx = piece_idx; idx < cumulative_max_parameters_.size(); ++idx) {
        if (idx == 0) {
            cumulative_max_parameters_[idx] = piece_ptrs_[idx]->max_parameter();
        } else {
            cumulative_max_parameters_[idx] =
                cumulative_max_parameters_[idx - 1] + piece_ptrs_[idx]->max_parameter();
        }
    }
}

template <typename T, unsigned int DIM>
void SingleParameterPiecewiseCurve<T, DIM>::scaleMaxParameters(T scaling_factor) {

    if (scaling_factor < 0) {
        throw std::invalid_argument("SingleParameterPiecewiseCurve::scaleMaxParameters: "
                                    "scaling_factor is negative. scaling_factor: " +
                                    std::to_string(scaling_factor));
    }

    for (std::size_t piece_idx = 0; piece_idx < numPieces(); ++piece_idx) {
        piece_ptrs_[piece_idx]->scaleMaxParameter(scaling_factor);
    }

    for (std::size_t piece_idx = 0; piece_idx < numPieces(); ++piece_idx) {
        if (piece_idx == 0) {
            cumulative_max_parameters_[piece_idx] = piece_ptrs_[piece_idx]->max_parameter();
        } else {
            cumulative_max_parameters_[piece_idx] =
                cumulative_max_parameters_[piece_idx - 1] + piece_ptrs_[piece_idx]->max_parameter();
        }
    }
}

template <typename T, unsigned int DIM>
T SingleParameterPiecewiseCurve<T, DIM>::max_parameter() const {
    if (cumulative_max_parameters_.empty()) {
        throw std::runtime_error("SingleParameterPiecewiseCurve::max_parameter: max_parameter is "
                                 "called when there are 0 pieces in the piecewise curve");
    }

    return cumulative_max_parameters_.back();
}

template <typename T, unsigned int DIM>
typename SingleParameterPiecewiseCurve<T, DIM>::VectorDIM
SingleParameterPiecewiseCurve<T, DIM>::eval(T parameter, uint64_t derivative_degree) const {
    if (piece_ptrs_.empty()) {
        throw std::runtime_error("SingleParameterPiecewiseCurve::eval: eval is called "
                                 "when there are 0 pieces in the piecewise curve");
    }

    const T max_param = this->max_parameter();

    if (parameter < 0 || parameter > max_param) {
        throw std::invalid_argument(
            "SingleParameterPiecewiseCurve::eval: parameter is out of range.");
    }

    typename std::vector<T>::const_iterator cumulative_parameter_iterator = std::lower_bound(
        cumulative_max_parameters_.begin(), cumulative_max_parameters_.end(), parameter);

    if (cumulative_parameter_iterator == cumulative_max_parameters_.end()) {
        throw std::runtime_error(
            "SingleParameterPiecewiseCurve::eval: std::lower_bound failed, but "
            "it should not have failed.");
    }

    std::size_t piece_idx =
        std::distance(cumulative_max_parameters_.begin(), cumulative_parameter_iterator);

    if (piece_idx == 0) {
        return piece_ptrs_[piece_idx]->eval(parameter, derivative_degree);
    } else {
        const T piece_max_parameter = piece_ptrs_[piece_idx]->max_parameter();
        return piece_ptrs_[piece_idx]->eval(
            std::min(piece_max_parameter, parameter - cumulative_max_parameters_[piece_idx - 1]),
            derivative_degree);
    }
}

template <typename T, unsigned int DIM>
typename SingleParameterPiecewiseCurve<T, DIM>::MaximumDerivativeMagnitude
SingleParameterPiecewiseCurve<T, DIM>::maximumDerivativeMagnitude(
    uint64_t derivative_degree) const {
    const T max_parameter = this->max_parameter();

    T maximum_derivative_magnitude = std::numeric_limits<T>::lowest();
    T maximum_derivative_magnitude_parameter = T(0.0);
    for (std::size_t piece_idx = 0; piece_idx < piece_ptrs_.size(); ++piece_idx) {

        const MaximumDerivativeMagnitude& magnitude_parameter_pair =
            piece_ptrs_[piece_idx]->maximumDerivativeMagnitude(derivative_degree);

        if (magnitude_parameter_pair.magnitude > maximum_derivative_magnitude) {
            maximum_derivative_magnitude = magnitude_parameter_pair.magnitude;
            maximum_derivative_magnitude_parameter =
                std::min(max_parameter,
                         (piece_idx == 0 ? T(0.0) : cumulative_max_parameters_.at(piece_idx - 1)) +
                             magnitude_parameter_pair.parameter);
        }
    }

    return MaximumDerivativeMagnitude{.magnitude = maximum_derivative_magnitude,
                                      .parameter = maximum_derivative_magnitude_parameter};
}

template class SingleParameterPiecewiseCurve<double, 3U>;

template class SingleParameterPiecewiseCurve<double, 2U>;

template class SingleParameterPiecewiseCurve<float, 3U>;

template class SingleParameterPiecewiseCurve<float, 2U>;

} // namespace splines