#pragma once

#include <cbf/detail/ConnectivityCBF.h>
#include <math/Geometry.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <mpc_cbf/controller/ConnectivityIMPCCBF.h>
#include <nlohmann/json.hpp>

namespace common {

/**
 * @brief Parse PiecewiseBezierParams from JSON configuration
 * @tparam T Numeric type (typically double)
 * @tparam DIM Dimension of the system
 * @param config_json JSON configuration object
 * @return PiecewiseBezierParams object
 */
template <typename T, unsigned int DIM>
mpc::PiecewiseBezierParams<T, DIM> parsePiecewiseBezierParams(const nlohmann::json& config_json) {
    size_t num_pieces = config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = config_json["bezier_params"]["num_control_points"];
    T piece_max_parameter = config_json["bezier_params"]["piece_max_parameter"];

    return {num_pieces, num_control_points, piece_max_parameter};
}

/**
 * @brief Parse MPCParams from JSON configuration
 * @tparam T Numeric type (typically double)
 * @param config_json JSON configuration object
 * @return MPCParams object
 */
template <typename T>
mpc::MPCParams<T> parseMPCParams(const nlohmann::json& config_json) {
    using Vector = math::Vector<T>;
    using VectorDIM = math::VectorDIM<T, 3>; // Assuming 3D for this implementation

    // MPC parameters
    T h = config_json["mpc_params"]["h"];
    T Ts = config_json["mpc_params"]["Ts"];
    int k_hor = config_json["mpc_params"]["k_hor"];
    T w_pos_err = config_json["mpc_params"]["mpc_tuning"]["w_pos_err"];
    T w_u_eff = config_json["mpc_params"]["mpc_tuning"]["w_u_eff"];
    int spd_f = config_json["mpc_params"]["mpc_tuning"]["spd_f"];

    // Physical limits
    Vector p_min = Vector::Zero(2);
    p_min << config_json["physical_limits"]["p_min"][0], config_json["physical_limits"]["p_min"][1];

    Vector p_max = Vector::Zero(2);
    p_max << config_json["physical_limits"]["p_max"][0], config_json["physical_limits"]["p_max"][1];

    VectorDIM v_min;
    v_min << config_json["physical_limits"]["v_min"][0], config_json["physical_limits"]["v_min"][1],
        config_json["physical_limits"]["v_min"][2];

    VectorDIM v_max;
    v_max << config_json["physical_limits"]["v_max"][0], config_json["physical_limits"]["v_max"][1],
        config_json["physical_limits"]["v_max"][2];

    VectorDIM a_min;
    a_min << config_json["physical_limits"]["a_min"][0], config_json["physical_limits"]["a_min"][1],
        config_json["physical_limits"]["a_min"][2];

    VectorDIM a_max;
    a_max << config_json["physical_limits"]["a_max"][0], config_json["physical_limits"]["a_max"][1],
        config_json["physical_limits"]["a_max"][2];

    return {h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, v_min, v_max, a_min, a_max}};
}

/**
 * @brief Parse IMPCParams from JSON configuration
 * @tparam T Numeric type (typically double)
 * @tparam DIM Dimension of the system
 * @param config_json JSON configuration object
 * @return IMPCParams object
 */
template <typename T, unsigned int DIM>
typename mpc_cbf::ConnectivityIMPCCBF<T, DIM>::IMPCParams
parseIMPCParams(const nlohmann::json& config_json) {
    // CBF parameters
    bool slack_mode = config_json["cbf_params"]["slack_mode"];
    T slack_cost = config_json["cbf_params"]["slack_cost"];
    T slack_decay_rate = config_json["cbf_params"]["slack_decay_rate"];
    int cbf_horizon = config_json["cbf_params"]["cbf_horizon"];
    int impc_iter = config_json["cbf_params"]["impc_iter"];
    
    return {cbf_horizon, impc_iter, slack_cost, slack_decay_rate, slack_mode};
}

/**
 * @brief Parse ConnectivityCBFParams from JSON configuration
 * @tparam T Numeric type (typically double)
 * @param config_json JSON configuration object
 * @return ConnectivityCBFParams object
 */
template <typename T>
cbf::ConnectivityCBFParams<T> parseConnectivityCBFParams(const nlohmann::json& config_json) {
    T d_min = config_json["cbf_params"]["d_min"];
    T d_max = config_json["cbf_params"]["d_max"];
    return {d_min, d_max};
}

/**
 * @brief Parse collision shape parameters and create AlignedBoxCollisionShape
 * @tparam T Numeric type (typically double)
 * @tparam DIM Dimension of the system
 * @param config_json JSON configuration object
 * @return Shared pointer to AlignedBoxCollisionShape
 */
template <typename T, unsigned int DIM>
std::shared_ptr<const math::AlignedBoxCollisionShape<T, DIM>>
parseCollisionShape(const nlohmann::json& config_json) {
    using VectorDIM = math::VectorDIM<T, DIM>;
    using AlignedBox = math::AlignedBox<T, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<T, DIM>;

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec << config_json["robot_params"]["collision_shape"]["aligned_box"][0],
        config_json["robot_params"]["collision_shape"]["aligned_box"][1],
        config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    return std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);
}

} // namespace common
