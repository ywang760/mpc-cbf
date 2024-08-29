//
// Created by lishuo on 1/31/24.
//

#ifndef SPLINES_SINGLEPARAMETERCURVE_H
#define SPLINES_SINGLEPARAMETERCURVE_H

#include <math/Types.h>
#include <memory>
#include <cstdint>

namespace splines {

    template<typename T, unsigned int DIM>
    class SingleParameterCurve {
    public:
        using VectorDIM = math::VectorDIM<T, DIM>;

        virtual ~SingleParameterCurve() = default;

        // get the maximum parameter of the curve. curve is assumed to be defined in
        // parameter range [0, max_parameter()]
        virtual T max_parameter() const = 0;

        // set the max parameter of the curve.
        virtual void set_max_parameter(T max_parameter) = 0;

        // scale max parameter by multiplying it with the scaling factor
        virtual void scaleMaxParameter(T scaling_factor) = 0;

        // Evaluate {derivative_degree}^{th} derivative of the curve at parameter
        virtual VectorDIM eval(T parameter, uint64_t derivative_degree) const = 0;

        /**
         * @brief Maximum derivative magnitude return value struct
         *
         */
        struct MaximumDerivativeMagnitude {
            // Maximum magnitude
            T magnitude;

            // Parameter at which the maximum occurs
            T parameter;
        };

        /**
         * @brief Returns the maximum magnitude of the derivative_degree^{th}
         * derivative of the curve as well as the parameter at which the maximum
         * occurs.
         *
         * @param derivative_degree The degree of the derivative.
         * @return MaximumDerivativeMagnitude The maximum magnitude
         * - parameter pair
         */
        virtual MaximumDerivativeMagnitude
        maximumDerivativeMagnitude(uint64_t derivative_degree) const = 0;

    };

} // splines

#endif //SPLINES_SINGLEPARAMETERCURVE_H
