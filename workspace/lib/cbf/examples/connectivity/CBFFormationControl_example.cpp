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

// Define state dimensions and types
constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;

int main(int argc, char *argv[])
{
    // Type aliases for better readability
    using ConnectivityCBF = cbf::ConnectivityCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using ConnectivityControl = cbf::ConnectivityControl<double, DIM>;
    using json = nlohmann::json;
    using Matrix = math::Matrix<double>;
    using Vector = math::Vector<double>;

    spdlog::set_pattern("[%H:%M:%S] [%^%l%$] [%s:%# %!] %v");

    const std::string DF_CFG = "/usr/src/mpc-cbf/workspace/experiments/config/examples/robots2_1.json";
    const std::string DF_OUT = "/usr/src/mpc-cbf/workspace/experiments/results/formation/states.json";

    // Parse command-line arguments
    cxxopts::Options options(
        "simulation",
        "connectivity simulation");

    options.add_options()("config_file", "Path to experiment configuration file",
                          cxxopts::value<std::string>()->default_value(DF_CFG))
                         ("write_filename", "Write output JSON to this file",
                          cxxopts::value<std::string>()->default_value(DF_OUT));

    auto option_parse = options.parse(argc, argv);
    // Load experiment configuration
    SPDLOG_INFO("Starting CBF Formation Control Example...");

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
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    int num_neighbors = experiment_config_json["tasks"]["so"].size() - 1;
    for (size_t i = 0; i < num_robots; ++i)
    {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], so_json[i][2];
        init_state.vel_ << VectorDIM::Zero();
        init_states.push_back(init_state);
        Vector current_state(6);
        current_state << init_state.pos_, init_state.vel_;
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
    }

    // init neighbor_ids
    std::vector<std::vector<size_t>> neighbor_ids(num_robots, std::vector<size_t>(num_robots - 1));
    for (size_t i = 0; i < num_robots; ++i)
    {
        size_t neighbor_idx = 0;
        for (size_t j = 0; j < num_robots; ++j)
        {
            if (i == j)
            {
                continue;
            }
            neighbor_ids[i][neighbor_idx] = j;
            neighbor_idx += 1;
        }
    }

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = Ts;
    states["Ts"] = Ts;

    // === Create one PID controller per robot ===
    std::vector<math::PID<double, 3U>> pid_controllers;
    math::PIDParams<double, 3U> pid_params{10.0, 0.1, 5.0, Ts}; // 你可以根据需要调参数

    for (size_t i = 0; i < num_robots; ++i) {
        pid_controllers.emplace_back(pid_params);
    }

    // Main simulation loop
    int loop_idx = 0;
    while (loop_idx < 100)
    {
        // Print out current loop_idx if loop_idx is a multiple of 10
        if (loop_idx % 10 == 0)
        {
            SPDLOG_INFO("Timestep: {}", loop_idx);
        }
        // Process each robot in the simulation
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx)
        {
            // These vectors will store estimated positions of other robots
            std::vector<VectorDIM> other_robot_positions;
            other_robot_positions.reserve(num_neighbors);

            // Process each neighbor robot
            for (int j = 0; j < num_robots - 1; ++j)
            {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                const VectorDIM &neighbor_pos = init_states.at(neighbor_id).pos_;
                VectorDIM extended_pos;
                extended_pos << neighbor_pos(0), neighbor_pos(1), 0;
                other_robot_positions.push_back(extended_pos);
            }

            // Compute desired control using spring control toward target

            // Uncomment the following 2 lines if you want to use critically damped spring control
            // const VectorDIM &target_pos = target_positions.at(robot_idx);
            // VectorDIM desired_u = math::criticallyDampedSpringControl<double, DIM>(init_states.at(robot_idx), target_pos, 0.5);

            VectorDIM ref_pos = target_positions.at(robot_idx);
            VectorDIM ref_vel = VectorDIM::Zero();  // 目标速度为0
            VectorDIM ref_acc = VectorDIM::Zero();  // 目标加速度为0
            VectorDIM desired_u = pid_controllers.at(robot_idx).control(
                init_states.at(robot_idx),
                ref_pos,
                ref_vel,
                ref_acc
            );
            
            // Apply CBF to modify control for safety and connectivity
            ConnectivityControl connectivity_control(connectivity_cbf, num_neighbors, slack_mode, slack_cost, slack_decay_rate);
            VectorDIM cbf_u;
            bool success = connectivity_control.optimize(cbf_u, desired_u, init_states.at(robot_idx),
                                                         other_robot_positions, a_min, a_max);
            if (!success)
            {
                SPDLOG_INFO("Optimization failed for robot {} at timestep {}", robot_idx, loop_idx);
                cbf_u = VectorDIM::Zero(); // Use zero control if optimization fails
            }
            //  else {
            //     auto cbf_ptr = connectivity_control.getCBF(); // 新增接口访问 cbf_
            //     const auto& state = init_states.at(robot_idx);
            //     Eigen::VectorXd x_self(6);
            //     x_self << state.pos_, state.vel_;  // 前 3 是位置，后 3 是速度
            //     // ✅ 调用连通性约束
            //     Eigen::VectorXd Ac = -1.0 * cbf_ptr->getConnConstraints(x_self, other_robot_positions);
            //     auto Bc = cbf_ptr->getConnBound(x_self, other_robot_positions);
            //     std::cout << "[CBF CHECK] robot " << robot_idx << ", timestep " << loop_idx << std::endl;
            //     std::cout << "  desired_u = " << desired_u.transpose() << std::endl;
            //     std::cout << "  actual_u  = " << cbf_u.transpose() << std::endl;
            //     std::cout << "  Ac * desired_u = " << Ac.dot(desired_u) << ", Bc = " << Bc << std::endl;
            //     std::cout << "  Ac * actual_u  = " << Ac.dot(cbf_u)    << ", Bc = " << Bc << std::endl;
            //     if (Ac.dot(desired_u) > Bc)
            //         std::cout << "  ➤ CBF modified the control input.\n";
            //     else
            //         std::cout << "  ➤ CBF did NOT need to modify the control input.\n";
            // }

            // Apply control to robot model to get next state and add noise
            State next_state = pred_model_ptr->applyInput(init_states.at(robot_idx), cbf_u);
            Vector x_t = (Vector(6) << next_state.pos_, next_state.vel_).finished();
            x_t = math::addRandomNoise(x_t, pos_std, vel_std);
            current_states.at(robot_idx) = x_t;

            // Log robot state, desried control, and cbf control for visualization/analysis
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
            states["robots"][std::to_string(robot_idx)]["desired_u"].push_back({desired_u[0], desired_u[1], desired_u[2]});
            states["robots"][std::to_string(robot_idx)]["cbf_u"].push_back({cbf_u[0], cbf_u[1], cbf_u[2]});
        }

        // Update each robot's state for the next iteration
        for (size_t robot_idx = 0; robot_idx < num_robots; ++robot_idx)
        {
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) = current_states.at(robot_idx)(2);
            init_states.at(robot_idx).vel_(0) = current_states.at(robot_idx)(3);
            init_states.at(robot_idx).vel_(1) = current_states.at(robot_idx)(4);
            init_states.at(robot_idx).vel_(2) = current_states.at(robot_idx)(5);
        }
        loop_idx += 1;
    }

    // Save simulation results to JSON file
    SPDLOG_INFO("Writing states to JSON file: {}", JSON_FILENAME);
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}