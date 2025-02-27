//
// Created by lishuo on 4/4/24.
//

#ifndef SPLINES_HELPERS_H
#define SPLINES_HELPERS_H

#include <nlohmann/json.hpp>
#include <fstream>
#include <splines/curves/SingleParameterPiecewiseCurve.h>
#include <splines/curves/Bezier.h>
#include <math/Types.h>

namespace splines {
    template <typename T, unsigned int DIM>
    nlohmann::json
    convertPiecewiseCurveToJson(
            splines::SingleParameterPiecewiseCurve<T, DIM>& piecewise_curve);

    template <typename T, unsigned int DIM>
    void writePiecewiseCurveToJson(
            splines::SingleParameterPiecewiseCurve<T, DIM>& piecewise_curve,
            std::string filename);

    template <typename T, unsigned int DIM>
    void writePiecewiseCurvesToJson(
            std::vector<splines::SingleParameterPiecewiseCurve<T, DIM>>& piecewise_curves,
            std::string filename);

} // splines

#endif //SPLINES_HELPERS_H
