//
// Created by lishuo on 8/17/24.
//

#include <iostream>
#include <memory>
#include <model/DoubleIntegratorXYYaw.h>
#include <nlohmann/json.hpp>
#include <fstream>

int main() {
    constexpr unsigned int DIM = 3U;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using StatePropagator = model::StatePropagator<double>;
    using State = model::State<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using json = nlohmann::json;

    double ts = 0.01;
    int K = 500;

    std::shared_ptr<DoubleIntegratorXYYaw> model_xyyaw_ptr = std::make_shared<DoubleIntegratorXYYaw>(ts);
    StatePropagator A0 = model_xyyaw_ptr->get_A0(K);
    StatePropagator Lambda = model_xyyaw_ptr->get_lambda(K);


    State initial_state;
    initial_state.pos_ << Vector::Zero(DIM); initial_state.vel_ << Vector::Zero(DIM);
    Vector x_0(6);
    x_0 << initial_state.pos_, initial_state.vel_;

    // define the control input
    Vector U = Vector::Zero(3*K);
    int idx = 0;
    for (size_t i = 0; i < K; ++i) {
        if (i < 100) {
            U[idx] = 1;
        }
        if (i >= 100 && i < 200) {
            U[idx] = -1;
            U[idx+1] = 1;
        }
        idx += 3;
    }

    std::cout << "Lambda pos matrix shape: (" << Lambda.pos_.rows() << ", " << Lambda.pos_.cols() << ")\n";
    std::cout << "Lambda vel matrix shape: (" << Lambda.vel_.rows() << ", " << Lambda.vel_.cols() << ")\n";
    std::cout << "A0 pos matrix shape: (" << A0.pos_.rows() << ", " << A0.pos_.cols() << ")\n";
    std::cout << "A0 vel matrix shape: (" << A0.vel_.rows() << ", " << A0.vel_.cols() << ")\n";

    Matrix x_pos_pred = A0.pos_ * x_0 + Lambda.pos_ * U;
    Matrix x_vel_pred = A0.vel_ * x_0 + Lambda.vel_ * U;

    // write to json
    std::string JSON_FILENAME = "XYYawStates.json";
    json states;
    states["dt"] = ts;
    int ts_idx = 0;
    for (size_t i = 0; i < K; ++i) {
        Vector x_t_pos = x_pos_pred.middleRows(ts_idx, 3);
        Vector x_t_vel = x_vel_pred.middleRows(ts_idx, 3);
        Vector x_t(6);
        x_t << x_t_pos, x_t_vel;
        states["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
        ts_idx += 3;
    }
    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}