//
// Created by lishuo on 8/24/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <mpc/controller/BezierMPC.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <nlohmann/json.hpp>

int main() {
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
    // load experiment config
    std::string experiment_config_filename = "../../../config/config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);

    size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
    double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];

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

    VectorDIM a_min;
    a_min << experiment_config_json["mpc_params"]["physical_limits"]["a_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][2];
    VectorDIM a_max;
    a_max << experiment_config_json["mpc_params"]["physical_limits"]["a_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][2];

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec <<
    experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
    experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
    experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
            std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);

    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
    MPCParams mpc_params = {h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, a_min, a_max}};

    std::string JSON_FILENAME = "../../../tools/XYYawStates.json";
    json states;
    states["dt"] = Ts;

    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h/Ts));
    StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h/Ts));
    // init MPC
    uint64_t bezier_continuity_upto_degree = 4;
    BezierMPCParams bezier_mpc_params = {piecewise_bezier_params, mpc_params};
//    BezierMPC bezier_mpc(bezier_mpc_params, pred_model_ptr, bezier_continuity_upto_degree);

    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
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

    SingleParameterPiecewiseCurve traj;
    int loop_idx = 0;
    while (loop_idx < 200) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            std::vector<VectorDIM> other_robot_positions;
            for (int j = 0; j < num_robots; ++j) {
                if (j==robot_idx) {
                    continue;
                }
                other_robot_positions.push_back(init_states.at(j).pos_);
            }

            Vector ref_positions(DIM*k_hor);
            // static target reference
            ref_positions = target_positions.at(robot_idx).replicate(k_hor, 1);

            BezierMPC bezier_mpc(bezier_mpc_params, pred_model_ptr, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr);
            bool success = bezier_mpc.optimize(traj, init_states.at(robot_idx), other_robot_positions, ref_positions);
            Vector U = bezier_mpc.generatorDerivativeControlInputs(2);
//            std::cout << "curve eval: " << traj.eval(1.5, 0) << "\n";

            // log down the optimized curve prediction
            double t = 0;
            while (t <= 1.5) {
                VectorDIM pred_pos = traj.eval(t, 0);
                states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx].push_back({pred_pos(0), pred_pos(1), pred_pos(2)});
                t += 0.05;
            }
            //


            Matrix x_pos_pred = exe_A0.pos_ * current_states.at(robot_idx) + exe_Lambda.pos_ * U;
            Matrix x_vel_pred = exe_A0.vel_ * current_states.at(robot_idx) + exe_Lambda.vel_ * U;
            assert(int(h/Ts)*DIM == x_pos_pred.rows());
            int ts_idx = 0;
            for (size_t i = 0; i < int(h / Ts); ++i) {
                Vector x_t_pos = x_pos_pred.middleRows(ts_idx, 3);
                Vector x_t_vel = x_vel_pred.middleRows(ts_idx, 3);
                Vector x_t(6);
                x_t << x_t_pos, x_t_vel;
                states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                ts_idx += 3;
                current_states.at(robot_idx) = x_t;
            }
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) = current_states.at(robot_idx)(2);
            init_states.at(robot_idx).vel_(0) = current_states.at(robot_idx)(3);
            init_states.at(robot_idx).vel_(1) = current_states.at(robot_idx)(4);
            init_states.at(robot_idx).vel_(2) = current_states.at(robot_idx)(5);
        }
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}

