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
    uint64_t bezier_continuity_upto_degree =
        config_json["bezier_params"]["bezier_continuity_upto_degree"];

    return {num_pieces, num_control_points, piece_max_parameter, bezier_continuity_upto_degree};
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

    // Validation: Critical MPC parameter relationships
    if (Ts > h) {
        throw std::invalid_argument("Control timestep Ts (" + std::to_string(Ts) +
                                    ") must be <= MPC timestep h (" + std::to_string(h) + ")");
    }
    if (h <= 0 || Ts <= 0) {
        throw std::invalid_argument("Time parameters h and Ts must be positive");
    }

    // Critical: h must be integer multiple of Ts for control discretization
    T ratio = h / Ts;
    if (std::abs(ratio - std::round(ratio)) > 1e-10) {
        throw std::invalid_argument("MPC timestep h (" + std::to_string(h) +
                                    ") must be an integer multiple of control timestep Ts (" +
                                    std::to_string(Ts) + ")");
    }
    if (spd_f > k_hor) {
        throw std::invalid_argument("Speed factor spd_f (" + std::to_string(spd_f) +
                                    ") must be <= prediction horizon k_hor (" +
                                    std::to_string(k_hor) + ")");
    }
    if (spd_f < 1) {
        throw std::invalid_argument("Speed factor spd_f must be at least 1");
    }
    if (k_hor < 1) {
        throw std::invalid_argument("Prediction horizon k_hor must be at least 1");
    }

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

    // Validation: IMPC parameter constraints
    if (cbf_horizon < 1) {
        throw std::invalid_argument("CBF horizon must be at least 1");
    }
    if (impc_iter < 1) {
        throw std::invalid_argument("IMPC iterations must be at least 1");
    }
    if (slack_mode && slack_cost <= 0) {
        throw std::invalid_argument("Slack cost must be positive when slack_mode is enabled");
    }
    if (slack_mode && (slack_decay_rate <= 0 || slack_decay_rate > 1)) {
        throw std::invalid_argument("Slack decay rate must be in (0,1] when slack_mode is enabled");
    }

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
    
    // Parse slack configuration - support both old and new format
    cbf::SlackConfig slack_config;
    
    if (config_json["cbf_params"].contains("slack_config")) {
        // New format with separate slack controls
        const auto& slack_json = config_json["cbf_params"]["slack_config"];
        slack_config.safety_slack = slack_json.value("safety_slack", false);
        slack_config.clf_slack = slack_json.value("clf_slack", false);
        slack_config.connectivity_slack = slack_json.value("connectivity_slack", false);
        slack_config.safety_slack_cost = slack_json.value("safety_slack_cost", 100000.0);
        slack_config.clf_slack_cost = slack_json.value("clf_slack_cost", 50000.0);
        slack_config.connectivity_slack_cost = slack_json.value("connectivity_slack_cost", 25000.0);
        slack_config.slack_decay_rate = slack_json.value("slack_decay_rate", 0.1);
        
        // Validate slack configuration
        if (slack_config.slack_decay_rate <= 0.0 || slack_config.slack_decay_rate > 1.0) {
            throw std::invalid_argument("slack_decay_rate must be in (0,1], got: " + std::to_string(slack_config.slack_decay_rate));
        }
        if (slack_config.safety_slack && slack_config.safety_slack_cost <= 0.0) {
            throw std::invalid_argument("safety_slack_cost must be positive when safety_slack is enabled");
        }
        if (slack_config.clf_slack && slack_config.clf_slack_cost <= 0.0) {
            throw std::invalid_argument("clf_slack_cost must be positive when clf_slack is enabled");
        }
        if (slack_config.connectivity_slack && slack_config.connectivity_slack_cost <= 0.0) {
            throw std::invalid_argument("connectivity_slack_cost must be positive when connectivity_slack is enabled");
        }
    } else if (config_json["cbf_params"].contains("slack_mode")) {
        // Backward compatibility with old format
        bool old_slack_mode = config_json["cbf_params"]["slack_mode"];
        if (old_slack_mode) {
            // Enable all slack types for backward compatibility
            slack_config.safety_slack = true;
            slack_config.clf_slack = true;
            slack_config.connectivity_slack = true;
            slack_config.safety_slack_cost = config_json["cbf_params"].value("slack_cost", 50000.0);
            slack_config.clf_slack_cost = config_json["cbf_params"].value("slack_cost", 50000.0);
            slack_config.connectivity_slack_cost = config_json["cbf_params"].value("slack_cost", 50000.0);
            slack_config.slack_decay_rate = config_json["cbf_params"].value("slack_decay_rate", 0.1);
        }
    }
    
    return {d_min, d_max, slack_config};
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

/**
 * @brief Validate cross-parameter relationships that span multiple parameter groups
 * @tparam T Numeric type (typically double)
 * @tparam DIM Dimension of the system
 * @param mpc_params MPC parameters object
 * @param bezier_params Bezier parameters object
 * @param impc_params IMPC parameters object
 * @throws std::invalid_argument if cross-parameter relationships are invalid
 */
template <typename T, unsigned int DIM>
inline void validateCrossParameterRelationships(
    const mpc::MPCParams<T>& mpc_params, const mpc::PiecewiseBezierParams<T, DIM>& bezier_params,
    const typename mpc_cbf::ConnectivityIMPCCBF<T, DIM>::IMPCParams& impc_params) {

    // Extract parameters from objects
    T h = mpc_params.h_;
    int k_hor = mpc_params.k_hor_;
    int cbf_horizon = impc_params.cbf_horizon_;

    size_t num_pieces = bezier_params.num_pieces_;
    T piece_max_parameter = bezier_params.piece_max_parameter_;

    // Critical: CBF horizon must be <= MPC prediction horizon
    if (cbf_horizon > k_hor) {
        throw std::invalid_argument("CBF horizon (" + std::to_string(cbf_horizon) +
                                    ") must be <= MPC prediction horizon k_hor (" +
                                    std::to_string(k_hor) + ")");
    }

    // Critical: Bezier curve parameter range must accommodate MPC sampling range
    T total_bezier_parameter = static_cast<T>(num_pieces) * piece_max_parameter;
    T max_mpc_parameter =
        static_cast<T>(k_hor - 1) * h; // From h_samples_ calculation in controller

    if (max_mpc_parameter > total_bezier_parameter) {
        throw std::invalid_argument(
            "MPC sampling range [0, " + std::to_string(max_mpc_parameter) +
            "] exceeds Bezier curve parameter range [0, " + std::to_string(total_bezier_parameter) +
            "]. Either reduce k_hor to <= " +
            std::to_string(static_cast<int>(total_bezier_parameter / h) + 1) +
            " or increase num_pieces/piece_max_parameter in bezier_params");
    }
}

} // namespace common
