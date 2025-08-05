//
// Created by lishuo on 8/24/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <mpc/controller/BezierMPC.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <nlohmann/json.hpp>
#include <cxxopts.hpp>
#include <fstream>
#include <spdlog/spdlog.h>
#include <common/logging.hpp>

int main(int argc, char *argv[])
{
    constexpr unsigned int DIM = 3U;
    using BezierMPC = mpc::BezierMPC<double, DIM>;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using State = model::State<double, DIM>;
    using StatePropagator = model::StatePropagator<double>;
    using json = nlohmann::json;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using BezierMPCParams = BezierMPC::Params;

    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;

    // Initialize logging with environment variable support
    auto logger = common::initializeLogging();

    const std::string DF_CFG = "/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/line.json";
    const std::string DF_OUT = "/usr/src/mpc-cbf/workspace/experiments/results/XYYawStates.json";

    // Parse command-line arguments
    cxxopts::Options options(
        "bezier_mpc_simulation",
        "Bezier MPC XY Yaw simulation");

    options.add_options()("config_file", "Path to experiment configuration file",
                          cxxopts::value<std::string>()->default_value(DF_CFG))
                         ("write_filename", "Write output JSON to this file",
                          cxxopts::value<std::string>()->default_value(DF_OUT))
                         ("max_steps", "Maximum number of simulation steps",
                          cxxopts::value<int>()->default_value("100"));
                          
    auto option_parse = options.parse(argc, argv);

    // Load experiment configuration
    logger->info("Starting Bezier MPC XYYaw Example...");

    std::string experiment_config_filename = option_parse["config_file"].as<std::string>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);

    size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
    double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];

    double h = experiment_config_json["mpc_params"]["h"];
    double Ts = experiment_config_json["mpc_params"]["h"]; // Use h for Ts as in migrated example
    int k_hor = experiment_config_json["mpc_params"]["k_hor"];
    double w_pos_err = experiment_config_json["mpc_params"]["mpc_tuning"]["w_pos_err"];
    double w_u_eff = experiment_config_json["mpc_params"]["mpc_tuning"]["w_u_eff"];
    int spd_f = experiment_config_json["mpc_params"]["mpc_tuning"]["spd_f"];
    Vector p_min = Vector::Zero(2);
    p_min << experiment_config_json["physical_limits"]["p_min"][0],
        experiment_config_json["physical_limits"]["p_min"][1];
    Vector p_max = Vector::Zero(2);
    p_max << experiment_config_json["physical_limits"]["p_max"][0],
        experiment_config_json["physical_limits"]["p_max"][1];

    VectorDIM a_min;
    a_min << experiment_config_json["physical_limits"]["a_min"][0],
        experiment_config_json["physical_limits"]["a_min"][1],
        experiment_config_json["physical_limits"]["a_min"][2];
    VectorDIM a_max;
    a_max << experiment_config_json["physical_limits"]["a_max"][0],
        experiment_config_json["physical_limits"]["a_max"][1],
        experiment_config_json["physical_limits"]["a_max"][2];

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec << experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
        std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);

    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
    MPCParams mpc_params = {h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, a_min, a_max}};

    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = Ts;
    states["Ts"] = Ts;

    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h / Ts));
    StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h / Ts));
    // init MPC
    uint64_t bezier_continuity_upto_degree = 4;
    BezierMPCParams bezier_mpc_params = {piecewise_bezier_params, mpc_params};

    // main loop
    // load the tasks
    std::vector<State> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i)
    {
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
    while (loop_idx < option_parse["max_steps"].as<int>())
    {
        // Print out current loop_idx if loop_idx is a multiple of 10
        if (loop_idx % 10 == 0)
        {
            logger->info("Timestep: {}", loop_idx);
        }
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx)
        {
            auto &current_state = current_states.at(robot_idx);

            // Collect positions of all other robots for collision avoidance
            std::vector<VectorDIM> other_robot_positions;
            for (int j = 0; j < num_robots; ++j)
            {
                if (j != robot_idx)
                {
                    other_robot_positions.push_back(current_states.at(j).pos_);
                }
            }

            // Create reference trajectory (static target for all horizon steps)
            Vector ref_positions = target_positions.at(robot_idx).replicate(k_hor, 1);

            // Solve MPC optimization problem
            BezierMPC bezier_mpc(bezier_mpc_params, pred_model_ptr, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr);
            bool success = bezier_mpc.optimize(traj, current_state, other_robot_positions, ref_positions);

            if (!success)
            {
                logger->warn("Optimization failed for robot {} at timestep {}", robot_idx, loop_idx);
            }
            // Extract control inputs from optimized trajectory
            Vector U = bezier_mpc.generatorDerivativeControlInputs(2);

            // Log predicted trajectory curve for visualization
            json pred_curve_data = json::array();
            for (double t = 0; t <= 1.5; t += 0.05)
            {
                VectorDIM pred_pos = traj.eval(t, 0);
                pred_curve_data.push_back({pred_pos(0), pred_pos(1), pred_pos(2)});
            }
            states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx] = std::move(pred_curve_data);

            // Apply control input to system dynamics
            Vector current_vec(6);
            current_vec << current_state.pos_, current_state.vel_;
            Matrix x_pos_pred = exe_A0.pos_ * current_vec + exe_Lambda.pos_ * U;
            Matrix x_vel_pred = exe_A0.vel_ * current_vec + exe_Lambda.vel_ * U;

            const int num_steps = int(h / Ts);
            assert(num_steps * DIM == x_pos_pred.rows());

            // Simulate forward and log all intermediate states
            Vector final_state;
            for (int i = 0; i < num_steps; ++i)
            {
                const int idx = i * DIM;
                Vector x_t_pos = x_pos_pred.middleRows(idx, DIM);
                Vector x_t_vel = x_vel_pred.middleRows(idx, DIM);

                // Log state at each timestep
                states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t_pos(0), x_t_pos(1), x_t_pos(2),
                                                                                 x_t_vel(0), x_t_vel(1), x_t_vel(2)});

                if (i == num_steps - 1)
                {
                    final_state.resize(2 * DIM);
                    final_state << x_t_pos, x_t_vel;
                }
            }

            // Update robot's current state for next iteration
            current_state.pos_ = final_state.head<DIM>();
            current_state.vel_ = final_state.tail<DIM>();
        }
        ++loop_idx;
    }

    // Save simulation results to JSON file
    logger->info("Writing states to JSON file: {}", JSON_FILENAME);
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}
