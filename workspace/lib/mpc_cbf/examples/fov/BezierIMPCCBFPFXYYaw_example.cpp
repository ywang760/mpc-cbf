//
// Created by lishuo on 9/21/24.
//

#include <math/Geometry.h>
#include <math/Random.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/FovBezierIMPCCBF.h>
#include <mpc_cbf/optimization/FovMPCCBFQPGenerator.h>

#include <cxxopts.hpp>
#include <fstream>
#include <nlohmann/json.hpp>

int main(int argc, char* argv[]) {
    constexpr unsigned int DIM = 3U;
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using FovBezierIMPCCBF = mpc_cbf::FovBezierIMPCCBF<double, DIM>;
    using State = model::State<double, DIM>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using FoVCBFParams = cbf::FoVCBFParams<double>;
    using BezierMPCCBFParams = mpc_cbf::FovMPCCBFQPOperations<double, DIM>::Params;
    using IMPCParams = mpc_cbf::FovBezierIMPCCBF<double, DIM>::IMPCParams;
    using IMPCCBFParams = mpc_cbf::FovBezierIMPCCBF<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    // cxxopt
    cxxopts::Options options("simulation", "fovmpc simulation");
    options.add_options()("instance_type", "instance type for simulations",
                          cxxopts::value<std::string>()->default_value("circle"))(
        "num_robots", "number of robots in the simulation",
        cxxopts::value<int>()->default_value(std::to_string(2)))(
        "fov", "fov degree", cxxopts::value<int>()->default_value(std::to_string(120)))(
        "slack_decay", "slack variable cost decay rate",
        cxxopts::value<double>()->default_value(std::to_string(0.1)))(
        "write_filename", "write to json filename",
        cxxopts::value<std::string>()->default_value(
            "/usr/src/mpc-cbf/workspace/experiments/results/states.json"));
    auto option_parse = options.parse(argc, argv);

    // load experiment config
    std::cout << "loading experiment settings...\n";
    std::string experiment_config_filename = "/usr/src/mpc-cbf/workspace/config/config.json";
    const int num_robots_parse = option_parse["num_robots"].as<int>();
    const int fov_beta_parse = option_parse["fov"].as<int>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    // piecewise bezier params
    size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
    double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];
    // mpc params
    double h = experiment_config_json["mpc_params"]["h"];
    double Ts = experiment_config_json["mpc_params"]["Ts"];
    int k_hor = experiment_config_json["mpc_params"]["k_hor"];
    double w_pos_err = experiment_config_json["mpc_params"]["mpc_tuning"]["w_pos_err"];
    double w_u_eff = experiment_config_json["mpc_params"]["mpc_tuning"]["w_u_eff"];
    int spd_f = experiment_config_json["mpc_params"]["mpc_tuning"]["spd_f"];
    Vector p_min = Vector::Zero(2);
    p_min << experiment_config_json["mpc_params"]["physical_limits"]["p_min"][0],
        experiment_config_json["mpc_params"]["physical_limits"]["p_min"][1];
    Vector p_max = Vector::Zero(2);
    p_max << experiment_config_json["mpc_params"]["physical_limits"]["p_max"][0],
        experiment_config_json["mpc_params"]["physical_limits"]["p_max"][1];
    // fov cbf params
    double fov_beta = double(fov_beta_parse) * M_PI / 180.0;
    std::cout << "fov_beta: " << double(fov_beta_parse) << "\n";
    double fov_Ds = experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

    // robot physical params
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

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec
        << experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
        experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
        std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);

    double pos_std = experiment_config_json["mpc_params"]["physical_limits"]["pos_std"];
    double vel_std = experiment_config_json["mpc_params"]["physical_limits"]["vel_std"];

    // create params
    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points,
                                                     piece_max_parameter};
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
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf =
        std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, v_min, v_max);
    // init bezier mpc-cbf
    uint64_t bezier_continuity_upto_degree = 3;
    int num_neighbors = experiment_config_json["tasks"]["so"].size() - 1;
    std::cout << "neighbor size: " << num_neighbors << "\n";
    BezierMPCCBFParams bezier_mpc_cbf_params = {piecewise_bezier_params, mpc_params,
                                                fov_cbf_params};
    int cbf_horizon = 2;
    int impc_iter = 2;
    double slack_cost = 1000;
    double slack_decay_rate = option_parse["slack_decay"].as<double>();
    std::cout << "slack_decay_rate: " << slack_decay_rate << "\n";
    bool slack_mode = true;
    IMPCParams impc_params = {cbf_horizon, impc_iter, slack_cost, slack_decay_rate, slack_mode};
    IMPCCBFParams impc_cbf_params = {bezier_mpc_cbf_params, impc_params};
    FovBezierIMPCCBF bezier_impc_cbf(impc_cbf_params, pred_model_ptr, fov_cbf,
                                     bezier_continuity_upto_degree, aligned_box_collision_shape_ptr,
                                     num_neighbors);

    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    assert(num_robots_parse == num_robots);
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
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

    std::vector<std::vector<size_t>> neighbor_ids(num_robots, std::vector<size_t>(num_robots - 1));
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

    // planning results
    std::vector<std::shared_ptr<SingleParameterPiecewiseCurve>> pred_traj_ptrs(num_robots);
    std::vector<double> traj_eval_ts(num_robots, 0);
    double pred_horizon = num_pieces * piece_max_parameter;

    double sim_runtime = 40;
    double sim_t = 0;
    int loop_idx = 0;
    while (sim_t < sim_runtime) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;
            for (int j = 0; j < num_robots - 1; ++j) {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                // for debug: fixed estimate
                other_robot_positions.push_back(init_states.at(neighbor_id).pos_);

                Matrix other_robot_cov(DIM, DIM);
                other_robot_cov << 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1;
                other_robot_covs.push_back(other_robot_cov);
                states["robots"][std::to_string(robot_idx)]["estimates_mean"]
                      [std::to_string(neighbor_id)]
                          .push_back({init_states.at(neighbor_id).pos_(0),
                                      init_states.at(neighbor_id).pos_(1),
                                      init_states.at(neighbor_id).pos_(2)});
                states["robots"][std::to_string(robot_idx)]["estimates_cov"]
                      [std::to_string(neighbor_id)]
                          .push_back({other_robot_cov(0), other_robot_cov(1), other_robot_cov(2),
                                      other_robot_cov(3), other_robot_cov(4), other_robot_cov(5),
                                      other_robot_cov(6), other_robot_cov(7), other_robot_cov(8)});
            }
            //            BezierIMPCCBF bezier_impc_cbf(bezier_impc_cbf_params, pred_model_ptr,
            //            fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr,
            //            impc_iter);
            bezier_impc_cbf.resetProblem();
            Vector ref_positions(DIM * k_hor);
            // static target reference
            VectorDIM converted_target_position = math::convertToClosestYaw<DIM>(
                init_states.at(robot_idx), target_positions.at(robot_idx));
            ref_positions = converted_target_position.replicate(k_hor, 1);
            std::vector<SingleParameterPiecewiseCurve> trajs;
            bool success =
                bezier_impc_cbf.optimize(trajs, init_states.at(robot_idx), other_robot_positions,
                                         other_robot_covs, ref_positions);
            if (!success) {
                std::cout << "QP failed at ts: " << sim_t;
                if (!trajs.empty()) {
                    std::cout << "; But initial planning is successful...";
                    pred_traj_ptrs.at(robot_idx) =
                        std::make_shared<SingleParameterPiecewiseCurve>(std::move(trajs.back()));
                    traj_eval_ts.at(robot_idx) = 0;
                }
                std::cout << "\n";
            }

            //           continuous control
            if (success) {
                pred_traj_ptrs.at(robot_idx) =
                    std::make_shared<SingleParameterPiecewiseCurve>(std::move(trajs.back()));
                traj_eval_ts.at(robot_idx) = 0;
            }

            // log down the prediction
            double t = 0;
            while (t <= pred_horizon) {
                VectorDIM pred_pos;
                double pred_t = traj_eval_ts.at(robot_idx) + t;
                if (pred_t > pred_traj_ptrs.at(robot_idx)->max_parameter()) {
                    pred_t = pred_traj_ptrs.at(robot_idx)->max_parameter();
                }
                pred_pos = pred_traj_ptrs.at(robot_idx)->eval(pred_t, 0);
                states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx][0].push_back(
                    {pred_pos(0), pred_pos(1), pred_pos(2)});
                t += 0.05;
            }

            // log down the trajectory
            double eval_t = 0;
            for (int t_idx = 1; t_idx <= int(h / Ts); ++t_idx) {
                eval_t = traj_eval_ts.at(robot_idx) + Ts * t_idx;
                if (eval_t > pred_traj_ptrs.at(robot_idx)->max_parameter()) {
                    eval_t = pred_traj_ptrs.at(robot_idx)->max_parameter();
                }

                Vector x_t_pos = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 0);
                x_t_pos(2) = math::convertYawInRange(x_t_pos(2));
                Vector x_t_vel = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 1);
                Vector x_t(6);
                x_t << x_t_pos, x_t_vel;
                x_t = math::addRandomNoise(x_t, pos_std, vel_std);
                states["robots"][std::to_string(robot_idx)]["states"].push_back(
                    {x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                current_states.at(robot_idx) = x_t;
            }
            traj_eval_ts.at(robot_idx) = eval_t;
        }
        // update the init_states
        for (size_t robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) =
                math::convertYawInRange(current_states.at(robot_idx)(2));
            init_states.at(robot_idx).vel_(0) = current_states.at(robot_idx)(3);
            init_states.at(robot_idx).vel_(1) = current_states.at(robot_idx)(4);
            init_states.at(robot_idx).vel_(2) = current_states.at(robot_idx)(5);
        }
        sim_t += h;
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    return 0;
}