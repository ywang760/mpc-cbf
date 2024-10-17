//
// Created by lishuo on 9/2/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <cbf/detail/cbf.h>
#include <cbf/controller/CBFControl.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>

constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;

math::VectorDIM<double, DIM> criticallyDampedSpringControl(const State& current_state, const VectorDIM& target, const double spring_constant) {
    VectorDIM Fs;
    Fs = spring_constant * target - current_state.pos_;
    VectorDIM Fd;
    Fd = -current_state.vel_ * 2 * sqrt(spring_constant);
    return Fs + Fd;
}

math::VectorDIM<double, DIM> rotateControlInputToBodyFrame(const VectorDIM& control_wf, double radian) {
    math::Matrix<double> R(3,3);
    R << cos(radian), sin(radian), 0,
            -sin(radian), cos(radian), 0,
            0, 0, 1;
    return R * control_wf;
}

math::VectorDIM<double, DIM> rotateControlInputToWorldFrame(const VectorDIM& control_bf, double radian) {
    math::Matrix<double> R(3,3);
    R << cos(radian), sin(radian), 0,
            -sin(radian), cos(radian), 0,
            0, 0, 1;
    return R.transpose() * control_bf;
}

int main() {
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using CBFControl = cbf::CBFControl<double, DIM>;
    using json = nlohmann::json;

    using Vector = math::Vector<double>;
    // load experiment config
    std::string experiment_config_filename = "../../../config/config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    double h = 0.1;

    double fov_beta = 120.0 * M_PI / 180.0;
    double fov_Ds = 3.0;
    double fov_Rs = 8.0;
    double vmax = 2.0;

    VectorDIM a_min;
    a_min << -5,-5,-5;
    VectorDIM a_max;
    a_max << 5,5,5;
    double VMAX = 1.;
    // json for record
    std::string JSON_FILENAME = "../../../tools/CBFXYYawStates.json";
    json states;
    states["dt"] = h;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, vmax);

    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < 1; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << -1,-1,0.785398;
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

    VectorDIM target_pos;
    target_pos << 4, 2.5, 0.0;
    // control loop
    int loop_idx = 0;
    while (loop_idx < 2000) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            Vector& current_state = current_states.at(robot_idx);
            Vector target(2);
            target << target_pos(0), target_pos(1);

            VectorDIM desired_u;
            desired_u << 0.2, 0.0, 0.0;
            // cbf control
            CBFControl cbf_control(fov_cbf);
            VectorDIM cbf_u;

            bool qp_success = cbf_control.optimize(cbf_u, desired_u, current_state, target, a_min, a_max);
            VectorDIM cbf_u_wf;
            if (!qp_success) {
                std::cout << "qp fail at: " << loop_idx << "\n";
            }

            State next_init_state = pred_model_ptr->applyInput(init_states.at(robot_idx), cbf_u);

            init_states.at(robot_idx) = next_init_state;

            current_states.at(robot_idx)(0) = init_states.at(robot_idx).pos_(0);
            current_states.at(robot_idx)(1) = init_states.at(robot_idx).pos_(1);
            current_states.at(robot_idx)(2) = init_states.at(robot_idx).pos_(2);
            current_states.at(robot_idx)(3) = init_states.at(robot_idx).vel_(0);
            current_states.at(robot_idx)(4) = init_states.at(robot_idx).vel_(1);
            current_states.at(robot_idx)(5) = init_states.at(robot_idx).vel_(2);

            // log the state
            Vector x_t(6);
            x_t = current_states.at(robot_idx);
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
        }
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}