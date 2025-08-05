//
// Created by lishuo on 2/13/24.
//

#ifndef SPLINES_SINGLEPARAMETERPIECEWISECURVE_H
#define SPLINES_SINGLEPARAMETERPIECEWISECURVE_H

#include <splines/curves/SingleParameterCurve.h>

namespace splines {

template <typename T, unsigned int DIM>
class SingleParameterPiecewiseCurve {

  public:
    using SingleParameterCurve = splines::SingleParameterCurve<T, DIM>;
    using VectorDIM = typename SingleParameterCurve::VectorDIM;
    using MaximumDerivativeMagnitude = typename SingleParameterCurve::MaximumDerivativeMagnitude;

    // adds piece to the piecewise curve. ownership of the piece_ptr is
    // transferred to *this
    void addPiece(std::unique_ptr<SingleParameterCurve>&& piece_ptr);

    // return the number of pieces in the piecewise curve
    std::size_t numPieces() const;

    // return a const reference to the given piece
    std::reference_wrapper<const SingleParameterCurve> getPiece(std::size_t piece_idx) const;

    // set the max parameter of the given piece
    void setMaxParameter(std::size_t piece_idx, T new_max_parameter);

    // scales the max parameter of all pieces by the scaling_factor
    void scaleMaxParameters(T scaling_factor);

    // get the maximum parameter of the curve. curve is assumed to be defined in
    // parameter range [0, max_parameter()].
    T max_parameter() const;

    // Evaluate {derivative_degree}^{th} derivative of the curve at parameter.
    // return status is not ok if there are 0 pieces in the piecewise curve or
    // parameter is out of range
    VectorDIM eval(T parameter, uint64_t derivative_degree) const;

    /**
         * @brief Returns the maximum magnitude of the derivative_degree^{th}
         * derivative of the curve as well as the parameter at which the maximum
         * occurs.
         *
         * @param derivative_degree The degree of the derivative.
         * @return MaximumDerivativeMagnitude The maximum magnitude
         * - parameter pair
         */
    MaximumDerivativeMagnitude maximumDerivativeMagnitude(uint64_t derivative_degree) const;

  private:
    // pointers to pieces
    std::vector<std::unique_ptr<SingleParameterCurve>> piece_ptrs_;

    // cumulative max parameters of pieces. cumulative_max_parameters_[i]
    // contains sum of max parameters of pieces [0, ..., i]
    std::vector<T> cumulative_max_parameters_;
};

} // namespace splines

#endif //SPLINES_SINGLEPARAMETERPIECEWISECURVE_H
