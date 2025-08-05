//
// Created by lishuo on 4/4/24.
//

#ifndef SPLINES_HELPERS_H
#define SPLINES_HELPERS_H

#include <fstream>
#include <math/Types.h>
#include <nlohmann/json.hpp>
#include <splines/curves/Bezier.h>
#include <splines/curves/SingleParameterPiecewiseCurve.h>

namespace splines {
template <typename T, unsigned int DIM>
nlohmann::json
convertPiecewiseCurveToJson(splines::SingleParameterPiecewiseCurve<T, DIM>& piecewise_curve);

template <typename T, unsigned int DIM>
void writePiecewiseCurveToJson(splines::SingleParameterPiecewiseCurve<T, DIM>& piecewise_curve,
                               std::string filename);

template <typename T, unsigned int DIM>
void writePiecewiseCurvesToJson(
    std::vector<splines::SingleParameterPiecewiseCurve<T, DIM>>& piecewise_curves,
    std::string filename);

} // namespace splines

#endif //SPLINES_HELPERS_H
