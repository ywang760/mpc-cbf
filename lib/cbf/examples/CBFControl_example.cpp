//
// Created by lishuo on 9/2/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <cbf/detail/cbf.h>
#include <cbf/controller/CBFControl.h>
#include <particle_filter/detail/particle_filter.h>
#include <nlohmann/json.hpp>
#include <cxxopts.hpp>
#include <fstream>
#include <cmath>

constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;

math::Vector<double> closestPointOnEllipse(const math::VectorDIM<double, 3U> &robot_pos,
                                           const math::Vector<double> &target_mean,
                                           const math::Matrix<double> &target_cov) {
    if (!isinf(target_cov(0, 0))) {
        Eigen::EigenSolver<math::Matrix<double>> es(target_cov.block(0, 0, 3U-1, 3U-1));
        math::Vector<double> eigenvalues = es.eigenvalues().real();
        math::Matrix<double> eigenvectors = es.eigenvectors().real();

        // s = 4.605 for 90% confidence interval
        // s = 5.991 for 95% confidence interval
        // s = 9.210 for 99% confidence interval
        double s = 4.605;
        double a = sqrt(s * eigenvalues(0)); // major axis
        double b = sqrt(s * eigenvalues(1)); // minor axis

        // a could be smaller than b, so swap them
        if (a < b)
        {
            double temp = a;
            a = b;
            b = temp;
        }

        int m = 0; // higher eigenvalue index
        int l = 1; // lower eigenvalue index
        if (eigenvalues(1) > eigenvalues(0))
        {
            m = 1;
            l = 0;
        }

        double theta = atan2(eigenvectors(1, m), eigenvectors(0, m)); // angle of the major axis wrt positive x-asis (ccw rotation)
        if (theta < 0.0) {
            theta += M_PI;
        } // angle in [0, 2pi]

        double slope = atan2(-target_mean(1) + robot_pos(1), -target_mean(0) + robot_pos(0));
        double x_n = target_mean(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
        double y_n = target_mean(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

        math::Vector<double> p_near(2);
        p_near << x_n, y_n;
        return p_near;
//        double dist = sqrt(pow(p_near(0) - robot_pos(0), 2) + pow(p_near(1) - robot_pos(1), 2));
//
//        if (isnan(dist)) {
//            dist = 5.0;
//            return dist;
//        }
//
//        // Check if robot is inside ellipse
//        double d = sqrt(pow(target_mean(0) - robot_pos(0), 2) + pow(target_mean(1) - robot_pos(1), 2));
//        double range = sqrt(pow(target_mean(0) - p_near(0), 2) + pow(target_mean(1) - p_near(1), 2));
//
//        if (d < range) {
//            return -dist;
//        } else {
//            return dist;
//        }
    } else {
        math::Vector<double> p_near(2);
        p_near << 0, 0;
        return p_near;
    }

//    return -5.0;
}

double convertYawInRange(double yaw) {
    assert(yaw < 2 * M_PI && yaw > -2 * M_PI);
    if (yaw > M_PI) {
        return yaw - 2 * M_PI;
    } else if (yaw < -M_PI) {
        return yaw + 2 * M_PI;
    } else {
        return yaw;
    }
}

math::VectorDIM<double, 3U> convertToClosestYaw(model::State<double, 3U>& state, const math::VectorDIM<double, 3U>& goal) {
    using VectorDIM = math::VectorDIM<double, 3U>;
    using Vector = math::Vector<double>;

//    std::cout << "start get the current_yaw\n";
    double current_yaw = state.pos_(2);
    // generate the candidate desire yaws
    Vector candidate_yaws(3);

//    std::cout << "start build candidate yaw\n";
//    std::cout << goal_(2) << "\n";
//    std::cout << M_PI << "\n";
    candidate_yaws << goal(2), goal(2) + 2 * M_PI, goal(2) - 2 * M_PI;

//    std::cout << "compute the offset\n";
    Vector candidate_yaws_offset(3);
    candidate_yaws_offset << std::abs(candidate_yaws(0) - current_yaw), std::abs(candidate_yaws(1) - current_yaw), std::abs(candidate_yaws(2) - current_yaw);


//    std::cout << "start to find the argmin\n";
    // Find the index of the minimum element
    int argmin_index = 0;
    double min_value = candidate_yaws_offset(0);

    for (int i = 1; i < candidate_yaws_offset.size(); ++i) {
        if (candidate_yaws_offset(i) < min_value) {
            min_value = candidate_yaws_offset(i);
            argmin_index = i;
        }
    }

    VectorDIM converted_goal;
    converted_goal << goal(0), goal(1), candidate_yaws(argmin_index);
//    std::cout << "converted_goal: " << converted_goal.transpose() << "\n";
    return converted_goal;

}

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

bool insideFOV(Eigen::VectorXd robot, Eigen::VectorXd target, double fov, double range)
{
    double yaw = robot(2);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
            -sin(yaw), cos(yaw), 0.0,
            0.0, 0.0, 1.0;
    Eigen::VectorXd t_local = R.block<2,2>(0,0) * (target.head(2) - robot.head(2));
    double dist = t_local.norm();
    double angle = abs(atan2(t_local(1), t_local(0)));
    if (angle <= 0.5*fov && dist <= range)
    {
        return true;
    } else
    {
        return false;
    }
}

int main(int argc, char* argv[]) {
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using CBFControl = cbf::CBFControl<double, DIM>;
    using ParticleFilter = pf::ParticleFilter;
    using json = nlohmann::json;

    using Matrix = math::Matrix<double>;
    using Vector = math::Vector<double>;

    // cxxopt
    cxxopts::Options options(
            "simulation",
            "fovmpc simulation");
    options.add_options()
            ("instance_type", "instance type for simulations",
             cxxopts::value<std::string>()->default_value("circle"))
            ("num_robots", "number of robots in the simulation",
             cxxopts::value<int>()->default_value(std::to_string(2)))
            ("fov", "fov degree",
             cxxopts::value<int>()->default_value(std::to_string(120)))
            ("write_filename", "write to json filename",
             cxxopts::value<std::string>()->default_value("../../../experiments/instances/results/circle2States.json"));
    auto option_parse = options.parse(argc, argv);

    // load experiment config
    std::cout << "loading experiment settings...\n";
    // std::string experiment_config_filename = "../../../config/config.json";
    std::string instance_type = option_parse["instance_type"].as<std::string>();
    std::string instance_path = instance_type+"_instances/";
    const int num_robots_parse = option_parse["num_robots"].as<int>();
    const int fov_beta_parse = option_parse["fov"].as<int>();
    std::string experiment_config_filename = "../../../experiments/instances/"+instance_path+instance_type+std::to_string(num_robots_parse)+"_config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
//    double h = experiment_config_json["mpc_params"]["h"];
    double Ts = experiment_config_json["mpc_params"]["h"];

    double fov_beta = double(experiment_config_json["fov_cbf_params"]["beta"]) * M_PI / 180.0;
    double fov_Ds = experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];
    VectorDIM v_min;
    v_min << experiment_config_json["mpc_params"]["physical_limits"]["v_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["v_min"][1],
            -0.2;
//    v_min << -20, -20, -20;
    VectorDIM v_max;
    v_max << experiment_config_json["mpc_params"]["physical_limits"]["v_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["v_max"][1],
            0.2;

//    v_max << 20, 20, 20;

    VectorDIM a_min;
    a_min << experiment_config_json["mpc_params"]["physical_limits"]["a_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][2];
//    a_min << -120, -120, -120;
    VectorDIM a_max;
    a_max << experiment_config_json["mpc_params"]["physical_limits"]["a_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][2];
//    a_max << 120, 120, 120;
//    double VMAX = 1.;

    // filter params
    int num_particles = 100;
    Matrix initCov = 1.0*Eigen::MatrixXd::Identity(DIM-1, DIM-1);
    Matrix processCov = 0.25*Eigen::MatrixXd::Identity(DIM-1, DIM-1);
    Matrix measCov = 0.05*Eigen::MatrixXd::Identity(DIM-1, DIM-1);

    // json for record
    std::string JSON_FILENAME = "../../../tools/CBFXYYawStates.json";
    json states;
    states["dt"] = Ts;
    states["Ts"] = Ts;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, v_min, v_max);

    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    assert(num_robots_parse==num_robots);
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], convertYawInRange(so_json[i][2]);
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
    // init filters
    std::vector<std::vector<ParticleFilter>> filters;
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

    for (size_t i = 0; i < num_robots; ++i){
        filters.emplace_back(num_robots-1);
        for (size_t j = 0; j < num_robots-1; ++j) {
            size_t neighbor_id = neighbor_ids.at(i).at(j);
            // init particle filter
            Vector neighbor_xy(DIM-1);
            neighbor_xy << init_states.at(neighbor_id).pos_(0), init_states.at(neighbor_id).pos_(1);
            filters[i][j].init(num_particles, neighbor_xy, initCov, processCov, measCov);
        }
        assert(filters[i].size() == num_robots-1);
    }
    assert(filters.size() == num_robots);

    // control loop
    int loop_idx = 0;
    while (loop_idx < 400) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;
            for (int j = 0; j < num_robots-1; ++j) {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                const VectorDIM& ego_pos = init_states.at(robot_idx).pos_;
                const VectorDIM& neighbor_pos = init_states.at(neighbor_id).pos_;
                ParticleFilter &filter = filters.at(robot_idx).at(j);
                // update filter estimate
                filter.predict();
                Vector weights = filter.getWeights();
                Matrix samples = filter.getParticles();
                for (int s = 0; s < num_particles; ++s) {
                    if (insideFOV(ego_pos, samples.col(s), fov_beta, fov_Rs)) {
                        weights[s] /= 10.0;
                    }
                }
                filter.setWeights(weights);

                // emulate vision
                if (insideFOV(ego_pos, neighbor_pos, fov_beta, fov_Rs)) {
                    Vector neighbor_xy(DIM-1);
                    neighbor_xy << neighbor_pos(0), neighbor_pos(1);
                    filter.update(neighbor_xy);
                }

                filter.resample();
                filter.estimateState();

                // update variables
                Vector estimate = filter.getState();
                VectorDIM extended_estimate;
                extended_estimate << estimate(0), estimate(1), 0;
                Matrix cov(DIM-1, DIM-1);
                cov = filter.getDistribution();
                other_robot_positions.push_back(extended_estimate);
                Matrix extended_cov(DIM, DIM);
                extended_cov << cov(0), cov(1), 0,
                        cov(2), cov(3), 0,
                        0,      0,      0;
                other_robot_covs.push_back(extended_cov);
                // log estimate
                states["robots"][std::to_string(robot_idx)]["estimates_mean"][std::to_string(neighbor_id)].push_back({extended_estimate(0), extended_estimate(1), extended_estimate(2)});
                states["robots"][std::to_string(robot_idx)]["estimates_cov"][std::to_string(neighbor_id)].push_back({extended_cov(0), extended_cov(1), extended_cov(2),
                                                                                                                     extended_cov(3), extended_cov(4), extended_cov(5),
                                                                                                                     extended_cov(6), extended_cov(7), extended_cov(8)});

                // for the closest point on ellipse visualization
                Vector p_near = closestPointOnEllipse(ego_pos, estimate, cov);
                states["robots"][std::to_string(robot_idx)]["p_near_ellipse"][std::to_string(neighbor_id)].push_back({p_near(0), p_near(1)});

//                // for debug: fixed estimate
//                other_robot_positions.push_back(init_states.at(neighbor_id).pos_);
//
//                Matrix other_robot_cov(DIM, DIM);
//                other_robot_cov << 0.1, 0, 0,
//                                   0, 0.1, 0,
//                                   0, 0, 0.1;
//                other_robot_covs.push_back(other_robot_cov);
            }

            // compute the desired control
//            const VectorDIM& target_pos = target_positions.at(robot_idx);
            const VectorDIM& target_pos = convertToClosestYaw(init_states.at(robot_idx), target_positions.at(robot_idx));
//            Vector& current_state = current_states.at(robot_idx);
//            Vector the_other_robot_2d_pos(2);
//            the_other_robot_2d_pos << the_other_robot_position(0), the_other_robot_position(1);

            VectorDIM desired_u = criticallyDampedSpringControl(init_states.at(robot_idx), target_pos, 1.);
            // cbf control
            CBFControl cbf_control(fov_cbf);
            VectorDIM cbf_u;
            std::cout << "robot id: " + std::to_string(robot_idx) +", state: " << init_states.at(robot_idx).pos_.transpose() << " " << init_states.at(robot_idx).vel_.transpose() << "\n";
            std::cout << "neighbor state: \n";
            for (size_t j = 0; j < other_robot_positions.size(); ++j) {
                std::cout << other_robot_positions.at(j).transpose() << "; ";
            }
            std::cout << "\n";
//            << init_states.at(robot_idx).pos_.transpose() << " " << init_states.at(robot_idx).vel_.transpose() << "\n";
            bool success = cbf_control.optimize(cbf_u, desired_u, init_states.at(robot_idx), other_robot_positions, a_min, a_max);
            if (!success) {
                std::cout << "optimization failed\n";
//                cbf_u = desired_u;
                cbf_u = VectorDIM::Zero();
            }
            std::cout << "optimized control: " << cbf_u.transpose() << "\n";
            State next_state = pred_model_ptr->applyInput(init_states.at(robot_idx), cbf_u);
            Vector x_t_pos = next_state.pos_;
            x_t_pos(2) = convertYawInRange(x_t_pos(2));
            Vector x_t_vel = next_state.vel_;
            Vector x_t(6);
            x_t << x_t_pos, x_t_vel;
            current_states.at(robot_idx) = x_t;

//            init_states.at(robot_idx) = next_init_state;
//
//            current_states.at(robot_idx)(0) = init_states.at(robot_idx).pos_(0);
//            current_states.at(robot_idx)(1) = init_states.at(robot_idx).pos_(1);
//            current_states.at(robot_idx)(2) = init_states.at(robot_idx).pos_(2);
//            current_states.at(robot_idx)(3) = init_states.at(robot_idx).vel_(0);
//            current_states.at(robot_idx)(4) = init_states.at(robot_idx).vel_(1);
//            current_states.at(robot_idx)(5) = init_states.at(robot_idx).vel_(2);
//
//            // log the state
//            Vector x_t(6);
//            x_t = current_states.at(robot_idx);
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
        }
        // update the init_states
        for (size_t robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) = convertYawInRange(current_states.at(robot_idx)(2));
//            init_states.at(robot_idx).pos_(2) = current_states.at(robot_idx)(2);
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