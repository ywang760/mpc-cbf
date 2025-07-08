//
// Created by Yichun on 6/27/2025
//

#include <model/DoubleIntegratorXYYaw.h>
#include <cbf/detail/ConnectivityCBF.h>
#include <cbf/controller/ConnectivityControl.h>
#include <nlohmann/json.hpp>
#include <cxxopts.hpp>
#include <fstream>
#include <cmath>
#include <math/Geometry.h>
#include <math/Controls.h>
#include <math/Random.h>
#include <spdlog/spdlog.h>
#include <common/logging.hpp>

int main(int argc, char *argv[])
{
    constexpr unsigned int DIM = 3U;
    // Type aliases for better readability
    using ConnectivityCBF = cbf::ConnectivityCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using ConnectivityControl = cbf::ConnectivityControl<double, DIM>;
    using json = nlohmann::json;
    using Matrix = math::Matrix<double>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using State = model::State<double, DIM>;

    // Initialize logging with environment variable support
    auto logger = common::initializeLogging();

    const std::string DF_CFG = "/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/line.json";
    const std::string DF_OUT = "/usr/src/mpc-cbf/workspace/experiments/results/states.json";

    // Parse command-line arguments
    cxxopts::Options options(
        "simulation",
        "connectivity simulation");

    options.add_options()("config_file", "Path to experiment configuration file",
                          cxxopts::value<std::string>()->default_value(DF_CFG))
                         ("write_filename", "Write output JSON to this file",
                          cxxopts::value<std::string>()->default_value(DF_OUT))
                         ("max_steps", "Maximum number of simulation steps",
                          cxxopts::value<int>()->default_value("100"));

    auto option_parse = options.parse(argc, argv);
    // Load experiment configuration
    logger->info("Starting CBF Formation Control Example...");

    // Configuration file path
    std::string experiment_config_filename = option_parse["config_file"].as<std::string>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);

    // Extract simulation parameters
    double Ts = experiment_config_json["mpc_params"]["h"]; // Time step

    // physical settings
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

    double pos_std = experiment_config_json["physical_limits"]["pos_std"];
    double vel_std = experiment_config_json["physical_limits"]["vel_std"];

    // init model
    // TODO: change the model to DoubleIntegratorXY (no yaw)
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);

    // init cbf
    double d_min = experiment_config_json["cbf_params"]["d_min"];
    double d_max = experiment_config_json["cbf_params"]["d_max"];
    std::shared_ptr<ConnectivityCBF> connectivity_cbf = std::make_shared<ConnectivityCBF>(d_min, d_max, v_min, v_max);

    // cbf controller setting
    bool slack_mode = experiment_config_json["cbf_params"]["slack_mode"];
    double slack_cost = experiment_config_json["cbf_params"]["slack_cost"];
    double slack_decay_rate = experiment_config_json["cbf_params"]["slack_decay_rate"];

    // load the tasks
    std::vector<State> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    size_t num_neighbors = num_robots - 1;
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

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = Ts;
    states["Ts"] = Ts;

    // === Create one PID controller per robot ===
    std::vector<math::PID<double, 3U>> pid_controllers;
    math::PIDParams<double, 3U> pid_params{
        experiment_config_json["pid_params"]["kp"],
        experiment_config_json["pid_params"]["ki"],
        experiment_config_json["pid_params"]["kd"],
        Ts};
        
    for (size_t i = 0; i < num_robots; ++i)
    {
        pid_controllers.emplace_back(pid_params);
    }

    // Main simulation loop
    int loop_idx = 0;
    while (loop_idx < option_parse["max_steps"].as<int>())
    {
        // Print out current loop_idx if loop_idx is a multiple of 10
        if (loop_idx % 10 == 0)
        {
            logger->info("Timestep: {}", loop_idx);
        }
        // Process each robot in the simulation
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx)
        {
            // Compute desired control using spring control toward target

            // Uncomment the following 2 lines if you want to use critically damped spring control
            const VectorDIM &target_pos = target_positions.at(robot_idx);
            VectorDIM desired_u = math::criticallyDampedSpringControl<double, DIM>(current_states.at(robot_idx), target_pos, 0.5);

            // VectorDIM ref_pos = target_positions.at(robot_idx);
            // VectorDIM ref_vel = VectorDIM::Zero();  // 目标速度为0
            // VectorDIM ref_acc = VectorDIM::Zero();  // 目标加速度为0
            // VectorDIM desired_u = pid_controllers.at(robot_idx).control(
            //     current_states.at(robot_idx),
            //     ref_pos,
            //     ref_vel,
            //     ref_acc
            // );

            // Apply CBF to modify control for safety and connectivity
            ConnectivityControl connectivity_control(connectivity_cbf, num_neighbors, slack_mode, slack_cost, slack_decay_rate);
            VectorDIM cbf_u;
            bool success = connectivity_control.optimize(cbf_u, desired_u, current_states, robot_idx, a_min, a_max);
            if (!success)
            {
                logger->warn("Optimization failed for robot {} at timestep {}", robot_idx, loop_idx);
                cbf_u = VectorDIM::Zero(); // Fallback to zero control if optimization fails
            }

            // Apply control to robot model to get next state and add noise
            State next_state = pred_model_ptr->applyInput(current_states.at(robot_idx), cbf_u);
            next_state = math::addRandomNoise<double, DIM>(next_state, pos_std, vel_std);
            current_states.at(robot_idx) = next_state;

            // Log robot state, desired control, and cbf control for visualization/analysis
            states["robots"][std::to_string(robot_idx)]["states"].push_back({next_state.pos_(0), next_state.pos_(1), next_state.pos_(2),
                                                                             next_state.vel_(0), next_state.vel_(1), next_state.vel_(2)});
            states["robots"][std::to_string(robot_idx)]["desired_u"].push_back({desired_u[0], desired_u[1], desired_u[2]});
            states["robots"][std::to_string(robot_idx)]["cbf_u"].push_back({cbf_u[0], cbf_u[1], cbf_u[2]});
        }

        loop_idx += 1;
    }
    
    // Save simulation results to JSON file
    logger->info("Writing states to JSON file: {}", JSON_FILENAME);
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    return 0;
}