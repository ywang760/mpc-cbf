//
// Created by lishuo on 9/21/24.
//

#include <cbf/detail/FovCBF.h>
#include <cxxopts.hpp>
#include <fstream>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/BezierMPCCBF.h>
#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

int main(int argc, char* argv[]) {
    constexpr unsigned int DIM = 3U;
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using BezierMPCCBF = mpc_cbf::BezierMPCCBF<double, DIM>;
    using State = model::State<double, DIM>;
    using StatePropagator = model::StatePropagator<double>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using FoVCBFParams = cbf::FoVCBFParams<double>;
    using BezierMPCCBFParams = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    // Initialize logging
    auto logger = spdlog::default_logger();

    const std::string DF_CFG = "/usr/src/mpc-cbf/workspace/config/config_new.json";
    const std::string DF_OUT =
        "/usr/src/mpc-cbf/workspace/experiments/results/BezierMPCCBFXYYawStates.json";

    // Parse command-line arguments
    cxxopts::Options options("simulation", "Bezier MPC-CBF XYYaw simulation");

    options.add_options()("config_file", "Path to experiment configuration file",
                          cxxopts::value<std::string>()->default_value(DF_CFG))(
        "write_filename", "Write output JSON to this file",
        cxxopts::value<std::string>()->default_value(DF_OUT))(
        "max_steps", "Maximum number of simulation steps",
        cxxopts::value<int>()->default_value("200"));

    auto option_parse = options.parse(argc, argv);

    // Load experiment configuration
    logger->info("Starting Bezier MPC-CBF XYYaw Example...");
    std::string experiment_config_filename = option_parse["config_file"].as<std::string>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    // piecewise bezier params
    size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
    double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];
    uint64_t bezier_continuity_upto_degree =
        experiment_config_json["bezier_params"]["bezier_continuity_upto_degree"];
    // mpc params
    double h = experiment_config_json["mpc_params"]["h"];
    double Ts = experiment_config_json["mpc_params"]["Ts"];
    int k_hor = experiment_config_json["mpc_params"]["k_hor"];
    double w_pos_err = experiment_config_json["mpc_params"]["mpc_tuning"]["w_pos_err"];
    double w_u_eff = experiment_config_json["mpc_params"]["mpc_tuning"]["w_u_eff"];
    int spd_f = experiment_config_json["mpc_params"]["mpc_tuning"]["spd_f"];
    // fov cbf params
    double fov_beta = double(experiment_config_json["fov_cbf_params"]["beta"]) * M_PI / 180.0;
    double fov_Ds = experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

    Vector p_min = Vector::Zero(2);
    p_min << experiment_config_json["physical_limits"]["p_min"][0],
        experiment_config_json["physical_limits"]["p_min"][1];
    Vector p_max = Vector::Zero(2);
    p_max << experiment_config_json["physical_limits"]["p_max"][0],
        experiment_config_json["physical_limits"]["p_max"][1];
    VectorDIM v_min;
    v_min << experiment_config_json["physical_limits"]["v_min"][0],
        experiment_config_json["physical_limits"]["v_min"][1],
        experiment_config_json["physical_limits"]["v_min"][2];
    VectorDIM v_max;
    v_max << experiment_config_json["physical_limits"]["v_max"][0],
        experiment_config_json["physical_limits"]["v_max"][1],
        experiment_config_json["physical_limits"]["v_max"][2];
    VectorDIM a_min;
    a_min << experiment_config_json["physical_limits"]["a_min"][0],
        experiment_config_json["physical_limits"]["a_min"][1],
        experiment_config_json["physical_limits"]["a_min"][2];
    VectorDIM a_max;
    a_max << experiment_config_json["physical_limits"]["a_max"][0],
        experiment_config_json["physical_limits"]["a_max"][1],
        experiment_config_json["physical_limits"]["a_max"][2];

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec
        << experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
        std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);
    // create params
    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points,
                                                     piece_max_parameter, bezier_continuity_upto_degree};
    MPCParams mpc_params = {
        h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, v_min, v_max, a_min, a_max}};
    FoVCBFParams fov_cbf_params = {fov_beta, fov_Ds, fov_Rs};

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = h;
    states["Ts"] = Ts;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr =
        std::make_shared<DoubleIntegratorXYYaw>(h);
    std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr =
        std::make_shared<DoubleIntegratorXYYaw>(Ts);
    StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h / Ts));
    StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h / Ts));
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf =
        std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, v_min, v_max);
    // init bezier mpc-cbf
    uint64_t bezier_continuity_upto_degree = 4;
    BezierMPCCBFParams bezier_mpc_cbf_params = {piecewise_bezier_params, mpc_params,
                                                fov_cbf_params};

    // main loop
    // load the tasks - simplified state management
    std::vector<State> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
        // load initial positions and set current state to initial state
        State current_state;
        current_state.pos_ << so_json[i][0], so_json[i][1], so_json[i][2];
        current_state.vel_ << VectorDIM::Zero();
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
    }

    SingleParameterPiecewiseCurve traj;
    int loop_idx = 0;
    while (loop_idx < option_parse["max_steps"].as<int>()) {
        // Print out current loop_idx if loop_idx is a multiple of 10
        if (loop_idx % 10 == 0) {
            logger->info("Timestep: {}", loop_idx);
        }
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            std::vector<VectorDIM> other_robot_positions;
            for (int j = 0; j < num_robots; ++j) {
                if (j == robot_idx) {
                    continue;
                }
                other_robot_positions.push_back(current_states.at(j).pos_);
            }
            BezierMPCCBF bezier_mpc_cbf(bezier_mpc_cbf_params, pred_model_ptr, fov_cbf,
                                        aligned_box_collision_shape_ptr);

            Vector ref_positions(DIM * k_hor);
            // static target reference
            ref_positions = target_positions.at(robot_idx).replicate(k_hor, 1);

            bool success = bezier_mpc_cbf.optimize(traj, current_states.at(robot_idx),
                                                   other_robot_positions, ref_positions);
            if (!success) {
                logger->warn("Optimization failed for robot {} at timestep {}", robot_idx,
                             loop_idx);
            }
            Vector U = bezier_mpc_cbf.generatorDerivativeControlInputs(2);

            // log down the optimized curve prediction
            double t = 0;
            while (t <= 1.5) {
                VectorDIM pred_pos = traj.eval(t, 0);
                states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx].push_back(
                    {pred_pos(0), pred_pos(1), pred_pos(2)});
                t += 0.05;
            }

            // Multi-rate discrete control: execute control at higher frequency Ts within MPC
            // timestep h Convert State to Vector format for matrix operations
            Vector current_state_vec(6);
            current_state_vec << current_states.at(robot_idx).pos_,
                current_states.at(robot_idx).vel_;

            Matrix x_pos_pred = exe_A0.pos_ * current_state_vec + exe_Lambda.pos_ * U;
            Matrix x_vel_pred = exe_A0.vel_ * current_state_vec + exe_Lambda.vel_ * U;
            assert(int(h / Ts) * DIM == x_pos_pred.rows());
            int ts_idx = 0;
            Vector final_state_vec(6);
            for (size_t i = 0; i < int(h / Ts); ++i) {
                Vector x_t_pos = x_pos_pred.middleRows(ts_idx, 3);
                Vector x_t_vel = x_vel_pred.middleRows(ts_idx, 3);
                Vector x_t(6);
                x_t << x_t_pos, x_t_vel;
                states["robots"][std::to_string(robot_idx)]["states"].push_back(
                    {x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                ts_idx += 3;
                final_state_vec = x_t;
            }
            // Update current state from the discrete control execution
            current_states.at(robot_idx).pos_(0) = final_state_vec(0);
            current_states.at(robot_idx).pos_(1) = final_state_vec(1);
            current_states.at(robot_idx).pos_(2) = final_state_vec(2);
            current_states.at(robot_idx).vel_(0) = final_state_vec(3);
            current_states.at(robot_idx).vel_(1) = final_state_vec(4);
            current_states.at(robot_idx).vel_(2) = final_state_vec(5);
        }
        loop_idx += 1;
    }

    // Save simulation results to JSON file
    logger->info("Writing states to JSON file: {}", JSON_FILENAME);
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    return 0;
}
