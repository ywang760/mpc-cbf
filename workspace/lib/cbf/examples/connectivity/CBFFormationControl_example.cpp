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
#include <PID3D.h>


// Define state dimensions and types
constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;

int main(int argc, char* argv[]) {
    // Type aliases for better readability
    using ConnectivityCBF = cbf::ConnectivityCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using ConnectivityControl = cbf::ConnectivityControl<double, DIM>;
    using json = nlohmann::json;
    using Matrix = math::Matrix<double>;
    using Vector = math::Vector<double>;

    // Parse command-line arguments
    cxxopts::Options options(
            "simulation",
            "connectivity simulation");
    options.add_options()
            ("instance_type", "instance type for simulations",
             cxxopts::value<std::string>()->default_value("circle"))
            ("num_robots", "number of robots in the simulation",
             cxxopts::value<int>()->default_value(std::to_string(3)))
            ("slack_decay", "slack variable cost decay rate",
             cxxopts::value<double>()->default_value(std::to_string(0.1)))
            ("config_file", "path to experiment configuration file",              // TODO: change the config_file
             cxxopts::value<std::string>()->default_value("/usr/src/mpc-cbf/workspace/experiments/config/Baseline/r3_line.json"))
            ("write_filename", "write to json filename",
             cxxopts::value<std::string>()->default_value("/usr/src/mpc-cbf/workspace/experiments/results/Baseline/linear/linear_3.json"));
    auto option_parse = options.parse(argc, argv);

    // Load experiment configuration
    std::cout << "Loading experiment settings...\n";
    std::string instance_type = option_parse["instance_type"].as<std::string>();
    std::string instance_path = instance_type+"_instances/";
    const int num_robots_parse = option_parse["num_robots"].as<int>();

    // Configuration file path
    std::string experiment_config_filename = option_parse["config_file"].as<std::string>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);

    // Extract simulation parameters
    double Ts = experiment_config_json["mpc_params"]["h"]; // Time step

    // Velocity and acceleration constraints
    VectorDIM v_min;
    v_min << experiment_config_json["mpc_params"]["physical_limits"]["v_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["v_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["v_min"][2];

    VectorDIM v_max;
    v_max << experiment_config_json["mpc_params"]["physical_limits"]["v_max"][0],
        experiment_config_json["mpc_params"]["physical_limits"]["v_max"][1],
        experiment_config_json["mpc_params"]["physical_limits"]["v_max"][2];

    VectorDIM a_min;
    a_min << experiment_config_json["mpc_params"]["physical_limits"]["a_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][2];

    VectorDIM a_max;
    a_max << experiment_config_json["mpc_params"]["physical_limits"]["a_max"][0],
        experiment_config_json["mpc_params"]["physical_limits"]["a_max"][1],
        experiment_config_json["mpc_params"]["physical_limits"]["a_max"][2];

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = Ts;
    states["Ts"] = Ts;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    // init cbf
    // TODO: make these parameters configurable
    double min_dist = 0.2;
    double max_dist = 2.0;
    std::shared_ptr<ConnectivityCBF> connectivity_cbf = std::make_shared<ConnectivityCBF>(min_dist, max_dist, v_min, v_max);
    // cbf controller setting
    bool slack_mode = false; // TODO: change this to be configurable
    double slack_cost = 1000;
    double slack_decay_rate = option_parse["slack_decay"].as<double>();
    if (slack_mode)
    {
        std::cout << "slack mode is enabled\n";
        std::cout << "slack_cost: " << slack_cost << "\n";
        std::cout << "slack_decay: " << slack_decay_rate << "\n";
    }
    else
    {
        std::cout << "slack mode is disabled\n";
    }
    // physical settings
    double pos_std = experiment_config_json["mpc_params"]["physical_limits"]["pos_std"];
    double vel_std = experiment_config_json["mpc_params"]["physical_limits"]["vel_std"];
    
    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    assert(num_robots_parse==num_robots);
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    int num_neighbors = experiment_config_json["tasks"]["so"].size() - 1;
    for (size_t i = 0; i < num_robots; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], math::convertYawInRange(so_json[i][2]);
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
    std::vector<std::vector<size_t>> neighbor_ids(num_robots, std::vector<size_t>(num_robots-1));
    for (size_t i = 0; i < num_robots; ++i) {
        size_t neighbor_idx = 0;
        for (size_t j = 0; j < num_robots; ++j) {
            if (i == j) {
                continue;
            }
            neighbor_ids[i][neighbor_idx] = j;
            neighbor_idx += 1;
        }
    }

    std::vector<int> failed_timesteps;
    std::vector<int> failed_robot_ids;
    std::vector<double> failed_h_vals;
    std::vector<Eigen::RowVectorXd> failed_Ac;
    std::vector<double> failed_Bc;
    std::vector<VectorDIM> failed_u_desired;
    std::vector<VectorDIM> failed_positions;

    // === Create one PID3D controller per robot ===
    std::vector<pid::PID3D<double>> pid_controllers;
    pid::PIDParams<double> pid_params{1.0, 0.05, 1.0, Ts};  // 你可以根据需要调参数

    for (size_t i = 0; i < num_robots; ++i) {
        pid_controllers.emplace_back(pid_params);
    }

    // Main simulation loop
    int loop_idx = 0;
    double h_this_timestep = 0.0;
    int qp_success_count = 0;
    std::ofstream clear_log("cbf_h_log.txt");
    clear_log.close();
    while (loop_idx < 3000) {
        // Print out current loop_idx if loop_idx is a multiple of 10
        if (loop_idx % 10 == 0)
        {
            std::cout << "Timestep: " << loop_idx << "\n";
        }
        // Process each robot in the simulation
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            // These vectors will store estimated positions and covariances of other robots
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;

            // Process each neighbor robot
            for (int j = 0; j < num_robots-1; ++j) {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                const VectorDIM& ego_pos = init_states.at(robot_idx).pos_;
                const VectorDIM &neighbor_pos = init_states.at(neighbor_id).pos_;

                // Fixed estimate for debugging
                auto estimate = init_states.at(neighbor_id).pos_;
                Matrix cov(DIM - 1, DIM - 1);
                cov << 0.1, 0,
                    0, 0.1;

                // Extend the estimate to include yaw dimension (set to 0) and store them
                VectorDIM extended_estimate;
                extended_estimate << estimate(0), estimate(1), 0; // Add yaw dimension (set to 0)
                other_robot_positions.push_back(extended_estimate);

                // Extend covariance to include yaw dimension
                Matrix extended_cov(DIM, DIM);
                extended_cov << cov(0), cov(1), 0,
                        cov(2), cov(3), 0,
                        0,      0,      0;
                other_robot_covs.push_back(extended_cov);

                // Record estimates for logging/visualization
                states["robots"][std::to_string(robot_idx)]["estimates_mean"][std::to_string(neighbor_id)].push_back({extended_estimate(0), extended_estimate(1), extended_estimate(2)});
                states["robots"][std::to_string(robot_idx)]["estimates_cov"][std::to_string(neighbor_id)].push_back({extended_cov(0), extended_cov(1), extended_cov(2),
                                                                                                                     extended_cov(3), extended_cov(4), extended_cov(5),
                                                                                                                     extended_cov(6), extended_cov(7), extended_cov(8)});

                // Calculate closest point on uncertainty ellipse for visualization
                Vector p_near = math::closestPointOnEllipse(ego_pos, estimate, cov);
                states["robots"][std::to_string(robot_idx)]["p_near_ellipse"][std::to_string(neighbor_id)].push_back({p_near(0), p_near(1)});
            }

            // Compute desired control using spring control toward target
            //const VectorDIM &target_pos = math::convertToClosestYaw<DIM>(init_states.at(robot_idx), target_positions.at(robot_idx));
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
            auto cbf_ptr = connectivity_control.getCBF();
            const auto& state = init_states.at(robot_idx);
            Eigen::VectorXd x_self(6);
            x_self << state.pos_, state.vel_;  // 前3是位置，后3是速度
            double h_val = 0.0;
            // 获取 Ac, Bc
            Eigen::RowVectorXd Ac = -1.0 * cbf_ptr->getConnConstraints(x_self, other_robot_positions, h_val);
            auto Bc = cbf_ptr->getConnBound(x_self, other_robot_positions);
            bool success = connectivity_control.optimize(cbf_u, desired_u, state, other_robot_positions, other_robot_covs, a_min, a_max);
            if (!success) {
                std::cout << "Optimization failed for robot " << robot_idx << " at timestep " << loop_idx << "\n";
                cbf_u = VectorDIM::Zero(); // fallback
                //cbf_u = desired_u;  // fallback: use PID control if CBF fails
                failed_timesteps.push_back(loop_idx);
                failed_robot_ids.push_back(robot_idx);
                failed_h_vals.push_back(h_val);  // 记录失败时的连通性值
                // 保存调试数据
                failed_Ac.push_back(Ac);                    // Ac: RowVectorXd
                failed_Bc.push_back(Bc);                    // Bc: double or Vector
                failed_u_desired.push_back(desired_u);      // u_nominal
                failed_positions.push_back(state.pos_);     // 机器人当前位置
            } else {
                qp_success_count++;
                // ✅ 成功时也输出 Ac 和 h
                std::cout << "[DEBUG] Success: timestep " << loop_idx << ", robot " << robot_idx << "\n";
                //std::cout << "[DEBUG] Ac.rows() = " << Ac.rows() << ", Ac.cols() = " << Ac.cols() << ", Ac.size() = " << Ac.size() << std::endl;
                std::cout << "[DEBUG] Ac = " << Ac << std::endl;
                std::cout << "[DEBUG] Bc = " << Bc << std::endl;
                std::cout << "[DEBUG] h_val = " << h_val << std::endl;
                // 只在 robot 0 上记录一次 h_val（可以保留）
                if (robot_idx == 0) {
                    h_this_timestep = h_val;
                }
            }
            
            // Apply control to robot model to get next state
            State next_state = pred_model_ptr->applyInput(init_states.at(robot_idx), cbf_u);
           //State next_state = pred_model_ptr->applyInput(init_states.at(robot_idx), desired_u);

            // Extract position and velocity, normalize yaw angle
            Vector x_t_pos = next_state.pos_;
            x_t_pos(2) = math::convertYawInRange(x_t_pos(2));
            Vector x_t_vel = next_state.vel_;

            // Combine into full state vector
            Vector x_t(6);
            x_t << x_t_pos, x_t_vel;

            // Add process noise to simulate real-world uncertainty
            x_t = math::addRandomNoise(x_t, pos_std, vel_std);
            current_states.at(robot_idx) = x_t;

            // Log robot state for visualization/analysis
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
            // Print out robot state, where x_t[0] is x, x_t[1] is y, x_t[2] is yaw, x_t[3] is vx, x_t[4] is vy, x_t[5] is yaw rate
            // std::cout << "Robot " << robot_idx << " state: "
            //           << "x: " << x_t[0] << ", "
            //           << "y: " << x_t[1] << ", "
            //           << "yaw: " << x_t[2] << ", "
            //           << "vx: " << x_t[3] << ", "
            //           << "vy: " << x_t[4] << ", "
            //           << "yaw rate: " << x_t[5] << "\n";
        }

        // Update each robot's state for the next iteration
        for (size_t robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) = math::convertYawInRange(current_states.at(robot_idx)(2));
            init_states.at(robot_idx).vel_(0) = current_states.at(robot_idx)(3);
            init_states.at(robot_idx).vel_(1) = current_states.at(robot_idx)(4);
            init_states.at(robot_idx).vel_(2) = current_states.at(robot_idx)(5);
        }
        // std::ofstream fout("cbf_h_log.txt", std::ios::app);
        // if (fout.is_open()) {
        //     fout << "timestep " << loop_idx << ", h = " << h_this_timestep << std::endl;
        //     fout.close();
        // } else {
        //     std::cerr << "[Error] Unable to open cbf_h_log.txt for writing!" << std::endl;
        // }
        loop_idx += 1;
    }

    std::cout << "\n[CBF FAILURE SUMMARY]\n";
    for (size_t i = 0; i < failed_timesteps.size(); ++i) {
        std::cout << "Timestep: " << failed_timesteps[i]
                  << ", Robot: " << failed_robot_ids[i]
                  << ", h = " << failed_h_vals[i] << "\n";
        std::cout << "  Position   : " << failed_positions[i].transpose() << "\n";
        std::cout << "  Desired u  : " << failed_u_desired[i].transpose() << "\n";
        std::cout << "  Ac (row)   : " << failed_Ac[i] << "\n";
        std::cout << "  Bc         : " << failed_Bc[i] << "\n";
    }
    

    // Save simulation results to JSON file
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    std::cout << "Total successful QP optimizations: " << qp_success_count << " out of " << num_robots * loop_idx << std::endl;

    return 0;
}