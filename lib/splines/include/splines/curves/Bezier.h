//
// Created by lishuo on 1/31/24.
//

#ifndef SPLINES_BEZIER_H
#define SPLINES_BEZIER_H

#include <math/Types.h>
#include <splines/curves/SingleParameterCurve.h>
#include <splines/detail/BezierOperations.h>

namespace splines {

    template <typename T, unsigned int DIM>
    class Bezier : public SingleParameterCurve<T, DIM> {
    public:
        using Base = SingleParameterCurve<T, DIM>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using MaximumDerivativeMagnitude = typename Base::MaximumDerivativeMagnitude;

        // construct a bezier curve with the given max parameter and initial control
        // points.
        Bezier(T max_parameter, const std::vector<VectorDIM>& control_points);

        // append a control point to bezier curve
        void appendControlPoint(const VectorDIM& control_point);

        // get number of control points of the bezier curve. degree of the bezier
        // curve is numControlPoints() - 1
        std::size_t numControlPoints() const;

        // degree of the bezier curve. return status is not ok if degree is not
        // defined, i.e. when the bezier curve has no control points
        uint64_t degree() const;

        // get the control point with the given index. return status is not OK if
        // the index is out of range.
        const VectorDIM getControlPoint(std::size_t control_point_idx) const;

        // get all the control points of the bezier curve as a list
        const std::vector<VectorDIM> getControlPoints() const;

        T max_parameter() const override;
        void set_max_parameter(T max_parameter) override;
        void scaleMaxParameter(T scaling_factor) override;

        // Evaluate {derivative_degree}^{th} derivative of the curve at parameter.
        // return status is not ok if parameter is out of range [0, max_parameter_]
        // or curve does not have any control points
        virtual VectorDIM eval(
                T parameter, uint64_t derivative_degree) const override;

        /**
         * @brief Returns the maximum magnitude of the derivative_degree^{th}
         * derivative of the curve as well as the parameter at which the maximum
         * occurs.
         *
         * @param derivative_degree The degree of the derivative.
         * @return absl::StatusOr<MaximumDerivativeMagnitude> The maximum magnitude
         * - parameter pair
         */
        virtual MaximumDerivativeMagnitude
        maximumDerivativeMagnitude(uint64_t derivative_degree) const override;

    private:
        // max parameter of the curve
        T max_parameter_;

        // control points of the bezier curve
        std::vector<VectorDIM> control_points_;
    };

} // splines

#endif //SPLINES_BEZIER_H