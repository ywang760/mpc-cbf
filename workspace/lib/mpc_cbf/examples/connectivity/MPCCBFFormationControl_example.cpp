//
// Created by Yutong on 8/1/2025
//

#include <cbf/detail/ConnectivityCBF.h>
#include <math/Geometry.h>
#include <math/Random.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/ConnectivityIMPCCBF.h>

#include <common/logging.hpp>
#include <common/parsing.hpp>
#include <cxxopts.hpp>
#include <fstream>
#include <nlohmann/json.hpp>

int main(int argc, char* argv[]) {
    constexpr unsigned int DIM = 3U;
    using ConnectivityCBF = cbf::ConnectivityCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using ConnectivityIMPCCBF = mpc_cbf::ConnectivityIMPCCBF<double, DIM>;
    using State = model::State<double, DIM>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using ConnectivityCBFParams = cbf::ConnectivityCBFParams<double>;
    using ConnectivityMPCCBFParams = mpc_cbf::ConnectivityMPCCBFQPOperations<double, DIM>::Params;
    using IMPCParams = mpc_cbf::ConnectivityIMPCCBF<double, DIM>::IMPCParams;
    using IMPCCBFParams = mpc_cbf::ConnectivityIMPCCBF<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    // Initialize logging with environment variable support
    auto logger = common::initializeLogging();

    const std::string DF_CFG =
        "/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/line.json";
    const std::string DF_OUT = "/usr/src/mpc-cbf/workspace/experiments/results/states.json";

    // Parse command-line arguments
    cxxopts::Options options("simulation", "connectivity simulation");

    options.add_options()("config_file", "Path to experiment configuration file",
                          cxxopts::value<std::string>()->default_value(DF_CFG))(
        "write_filename", "Write output JSON to this file",
        cxxopts::value<std::string>()->default_value(DF_OUT))(
        "sim_runtime", "Simulation runtime in seconds",
        cxxopts::value<double>()->default_value("40"));
    auto option_parse = options.parse(argc, argv);

    // Load experiment configuration
    logger->info("Starting MPC CBF Formation Control Example...");

    // Configuration file path
    std::string experiment_config_filename = option_parse["config_file"].as<std::string>();
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);

    // Parse parameters using helper functions
    PiecewiseBezierParams piecewise_bezier_params =
        common::parsePiecewiseBezierParams<double, DIM>(experiment_config_json);
    MPCParams mpc_params = common::parseMPCParams<double>(experiment_config_json);
    IMPCParams impc_params = common::parseIMPCParams<double, DIM>(experiment_config_json);
    ConnectivityCBFParams connectivity_cbf_params =
        common::parseConnectivityCBFParams<double>(experiment_config_json);
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
        common::parseCollisionShape<double, DIM>(experiment_config_json);

    // Assemble connectivity cbf params
    ConnectivityMPCCBFParams connectivity_mpc_cbf_params = {piecewise_bezier_params, mpc_params,
                                                            connectivity_cbf_params};
    IMPCCBFParams impc_cbf_params = {connectivity_mpc_cbf_params, impc_params};

    // additional parameters
    double pos_std = experiment_config_json["physical_limits"]["pos_std"];
    double vel_std = experiment_config_json["physical_limits"]["vel_std"];

    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr =
        std::make_shared<DoubleIntegratorXYYaw>(mpc_params.h_);

    // init connectivity cbf
    std::shared_ptr<ConnectivityCBF> connectivity_cbf = std::make_shared<ConnectivityCBF>(
        connectivity_cbf_params.dmin_, connectivity_cbf_params.dmax_, mpc_params.limits_.v_min_,
        mpc_params.limits_.v_max_);

    // init connectivity mpc-cbf
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    size_t num_neighbors = num_robots - 1;
    ConnectivityIMPCCBF connectivity_impc_cbf(impc_cbf_params, pred_model_ptr, connectivity_cbf,
                                              aligned_box_collision_shape_ptr, num_neighbors);

    // load the tasks
    std::vector<State> current_states;
    std::vector<VectorDIM> target_positions;
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
        // load initial states directly into current_states
        State current_state;
        current_state.pos_ << so_json[i][0], so_json[i][1], so_json[i][2];
        current_state.vel_ << VectorDIM::Zero();
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
    }

    // planning results
    std::vector<std::shared_ptr<SingleParameterPiecewiseCurve>> pred_traj_ptrs(num_robots);
    std::vector<double> traj_eval_ts(num_robots, 0);
    double pred_horizon =
        piecewise_bezier_params.num_pieces_ * piecewise_bezier_params.piece_max_parameter_;

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = mpc_params.h_;
    states["Ts"] = mpc_params.Ts_;

    // main simulation loop
    double sim_runtime = option_parse["sim_runtime"].as<double>();
    double sim_t = 0;
    int loop_idx = 0;
    while (sim_t < sim_runtime) {
        if (loop_idx % 50 == 0) {
            logger->info("loop_idx: {}, sim_t: {:.3f} seconds", loop_idx, sim_t);
        }

        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            connectivity_impc_cbf.resetProblem();

            Vector ref_positions(DIM * mpc_params.k_hor_);
            ref_positions = target_positions.at(robot_idx).replicate(mpc_params.k_hor_, 1);

            std::vector<SingleParameterPiecewiseCurve> trajs;
            bool success = connectivity_impc_cbf.optimize(trajs, /*current_states=*/current_states,
                                                          /*self_idx=*/robot_idx, ref_positions);

            if (!success) {
                logger->warn("Optimization failed for robot {} at sim_t {:.3f}", robot_idx, sim_t);
                if (trajs.empty()) {
                    logger->error("No previous trajectory available for robot {}", robot_idx);
                } else {
                    logger->warn("Returning the last successful trajectory for robot {}",
                                 robot_idx);
                }
            }

            if (!trajs.empty()) {
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
            for (int t_idx = 1; t_idx <= int(mpc_params.h_ / mpc_params.Ts_); ++t_idx) {
                eval_t = traj_eval_ts.at(robot_idx) + mpc_params.Ts_ * t_idx;
                if (eval_t > pred_traj_ptrs.at(robot_idx)->max_parameter()) {
                    eval_t = pred_traj_ptrs.at(robot_idx)->max_parameter();
                }

                Vector x_t_pos = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 0);
                Vector x_t_vel = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 1);
                State next_state = {x_t_pos, x_t_vel};
                next_state = math::addRandomNoise<double, DIM>(next_state, pos_std, vel_std);
                current_states.at(robot_idx) = next_state;

                states["robots"][std::to_string(robot_idx)]["states"].push_back(
                    {next_state.pos_(0), next_state.pos_(1), next_state.pos_(2), next_state.vel_(0),
                     next_state.vel_(1), next_state.vel_(2)});
            }
            traj_eval_ts.at(robot_idx) = eval_t;
        }
        sim_t += mpc_params.h_;
        loop_idx += 1;
    }

    // write states to json
    logger->info("Writing states to JSON file: {}", JSON_FILENAME);
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    return 0;
}