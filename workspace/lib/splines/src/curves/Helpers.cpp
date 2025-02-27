//
// Created by lishuo on 4/4/24.
//

#include <splines/curves/Helpers.h>

namespace splines {
    template <typename T, unsigned int DIM>
    nlohmann::json splines::convertPiecewiseCurveToJson(splines::SingleParameterPiecewiseCurve<T, DIM> &piecewise_curve) {
        using json = nlohmann::json;
        using SingleParameterCurve = splines::SingleParameterCurve<T, DIM>;
        using Bezier = splines::Bezier<T, DIM>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        // json for all control_pts for all pieces
        json bezier_piecewise_json;

        // go through all pieces and get all the control points from that piece
        for (size_t piece_idx = 0; piece_idx < piecewise_curve.numPieces(); ++piece_idx) {
            const SingleParameterCurve &piece = piecewise_curve.getPiece(piece_idx);
            const auto &bezier_piece = dynamic_cast<const Bezier &>(piece);
            // parameter for one piece
            T parameter = bezier_piece.max_parameter();
            // json for all control points for one piece
            json piece_control_pts_json;
            for (size_t control_pt_index = 0; control_pt_index < bezier_piece.numControlPoints(); ++control_pt_index) {
                // json for one control point for one piece
                json control_pt_json;
                VectorDIM control_pt = bezier_piece.getControlPoint(control_pt_index);
                for (size_t d = 0; d < DIM; ++d) {
                    control_pt_json.push_back(control_pt(d));
                }
                piece_control_pts_json.push_back(control_pt_json);
            }
            bezier_piecewise_json["control_points"].push_back(piece_control_pts_json);
            bezier_piecewise_json["parameters"].push_back(parameter);
        }
        return bezier_piecewise_json;
    }

    template <typename T, unsigned int DIM>
    void splines::writePiecewiseCurveToJson(splines::SingleParameterPiecewiseCurve<T, DIM>& piecewise_curve,
                                        std::string filename) {
        using json = nlohmann::json;
        json bezier_piecewise_json = convertPiecewiseCurveToJson<T, DIM>(piecewise_curve);
        // write to json
        std::ofstream o(filename);
        o << std::setw(4) << bezier_piecewise_json << std::endl;
    }

    template <typename T, unsigned int DIM>
    void splines::writePiecewiseCurvesToJson(std::vector<splines::SingleParameterPiecewiseCurve<T, DIM>>& piecewise_curves,
                                            std::string filename) {
        using json = nlohmann::json;
        json all_bezier_piecewise_json;
        for (size_t i = 0; i < piecewise_curves.size(); ++i) {
            json bezier_piecewise_json = convertPiecewiseCurveToJson<T, DIM>(piecewise_curves.at(i));
            all_bezier_piecewise_json.push_back(bezier_piecewise_json);
        }
        // write to json
        std::ofstream o(filename);
        o << std::setw(4) << all_bezier_piecewise_json << std::endl;
    }


    template nlohmann::json splines::convertPiecewiseCurveToJson<double, 3U>(
            splines::SingleParameterPiecewiseCurve<double, 3U> &piecewise_curve);
    template nlohmann::json splines::convertPiecewiseCurveToJson<float, 3U>(
            splines::SingleParameterPiecewiseCurve<float, 3U> &piecewise_curve);
    template nlohmann::json splines::convertPiecewiseCurveToJson<double, 2U>(
            splines::SingleParameterPiecewiseCurve<double, 2U> &piecewise_curve);
    template nlohmann::json splines::convertPiecewiseCurveToJson<float, 2U>(
            splines::SingleParameterPiecewiseCurve<float, 2U> &piecewise_curve);


    template void splines::writePiecewiseCurveToJson<double, 3U>(
            splines::SingleParameterPiecewiseCurve<double, 3U>& piecewise_curve,
            std::string filename);
    template void splines::writePiecewiseCurveToJson<float, 3U>(
            splines::SingleParameterPiecewiseCurve<float, 3U>& piecewise_curve,
            std::string filename);
    template void splines::writePiecewiseCurveToJson<double, 2U>(
            splines::SingleParameterPiecewiseCurve<double, 2U>& piecewise_curve,
            std::string filename);
    template void splines::writePiecewiseCurveToJson<float, 2U>(
            splines::SingleParameterPiecewiseCurve<float, 2U>& piecewise_curve,
            std::string filename);

    template void splines::writePiecewiseCurvesToJson<double, 3U>(
            std::vector<splines::SingleParameterPiecewiseCurve<double, 3U>>& piecewise_curves,
            std::string filename);
    template void splines::writePiecewiseCurvesToJson<float, 3U>(
            std::vector<splines::SingleParameterPiecewiseCurve<float, 3U>>& piecewise_curves,
            std::string filename);
    template void splines::writePiecewiseCurvesToJson<double, 2U>(
            std::vector<splines::SingleParameterPiecewiseCurve<double, 2U>>& piecewise_curves,
            std::string filename);
    template void splines::writePiecewiseCurvesToJson<float, 2U>(
            std::vector<splines::SingleParameterPiecewiseCurve<float, 2U>>& piecewise_curves,
            std::string filename);

} // splines