//
// Created by lishuo on 10/21/24.
//

#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/BezierIMPCCBF.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <nlohmann/json.hpp>

#include <ros/ros.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Pose.h>
#include <mavros_msgs/PositionTarget.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <tf2/utils.h>
#include <tf2/LinearMath/Quaternion.h>
#include <signal.h>

using std::placeholders::_1;
sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed
double ctrl_freq = 20.0;
geometry_msgs::PoseStamped current_pose;
mavros_msgs::State current_state;
double takeoff_time = 15;
double mission_time = 40;
double land_time = 5;

constexpr unsigned int DIM = 3U;

void state_cb(const mavros_msgs::State::ConstPtr& msg)
{
    current_state = *msg;
}


void pose_cb(const geometry_msgs::PoseStamped::ConstPtr& msg)
{
    current_pose = *msg;
}

class ControlNode
{
public:
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using BezierIMPCCBF = mpc_cbf::BezierIMPCCBF<double, DIM>;
    using State = model::State<double, DIM>;
    using StatePropagator = model::StatePropagator<double>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using FoVCBFParams = cbf::FoVCBFParams<double>;
    using BezierMPCCBFParams = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<double, DIM>::Params;
    using IMPCParams = mpc_cbf::BezierIMPCCBF<double, DIM>::IMPCParams;
    using IMPCCBFParams = mpc_cbf::BezierIMPCCBF<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    ControlNode() : nh_priv_("~")
    {
        // ---------------- ROS params -----------------
        this->nh_priv_.getParam("ROBOT_ID", ROBOT_ID);
        this->nh_priv_.getParam("NUM_TARGETS", NUM_TARGETS);
        this->nh_priv_.getParam("rate", rate);
        this->nh_priv_.getParam("CONFIG_FILENAME", CONFIG_FILENAME);


        // ---------- Subs and Pubs -------------------
        std::vector<size_t> neighbor_ids;
        for (size_t i = 0; i < NUM_TARGETS+1; ++i) {
            if (i == ROBOT_ID) {
                continue;
            }
            neighbor_ids.push_back(i);
        }
        assert(neighbor_ids.size() == NUM_TARGETS);


        mavros_state_sub_ = nh_.subscribe<mavros_msgs::State>("mavros/state", 10, state_cb);
        pose_sub_ = nh_.subscribe<geometry_msgs::PoseStamped>("mavros/local_position/pose", 10, pose_cb);
        arming_client_ = nh_.serviceClient<mavros_msgs::CommandBool>("mavros/cmd/arming");
        set_mode_client_ = nh_.serviceClient<mavros_msgs::SetMode>("mavros/set_mode");
        local_pos_pub_ = nh_.advertise<geometry_msgs::PoseStamped>("mavros/setpoint_position/local", 10);

        // Init params
        target_states_.resize(NUM_TARGETS);
        target_covs_.resize(NUM_TARGETS);
        target_states_[0] = VectorDIM::Zero();
        std::cout << "NUM_TARGETS: " << NUM_TARGETS << "\n";
//        std::cout << target_states_.size() << "\n";
//        state_.pos_ = VectorDIM::Zero();
        state_.vel_ = VectorDIM::Zero();
        // load experiment config
        std::string experiment_config_filename = CONFIG_FILENAME;
        std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
        json experiment_config_json = json::parse(experiment_config_fc);
        std::cout << "successfully load the json..." << "\n";

        // piecewise bezier params
        size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
        size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
        double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];
        // mpc params
        h_ = experiment_config_json["mpc_params"]["h"];
        Ts_ = experiment_config_json["mpc_params"]["Ts"];
        k_hor_ = experiment_config_json["mpc_params"]["k_hor"];
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
        double fov_beta = double(experiment_config_json["fov_cbf_params"]["beta"]) * M_PI / 180.0;
        double fov_Ds = experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0];
        double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

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
        aligned_box_collision_vec <<
                                  experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
                experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
                experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
        AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
        std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
                std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);

        // take off state
        takeoff_pose_.header.frame_id = "map";
        takeoff_pose_.pose.position.x = experiment_config_json["tasks"]["so"][ROBOT_ID][0];
        takeoff_pose_.pose.position.y = experiment_config_json["tasks"]["so"][ROBOT_ID][1];
        takeoff_pose_.pose.position.z = z_;
        double yaw = experiment_config_json["tasks"]["so"][ROBOT_ID][2];
//        std::cout << "init yaw: " << yaw << "\n";
        tf2::Quaternion q;
        q.setEuler(0, 0, yaw);
        takeoff_pose_.pose.orientation = tf2::toMsg(q);
//        std::cout << "init quaternion: " << q.x() << "," << q.y() << "," << q.z() << "," << q.w() << "\n";

        // create params
        PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
        MPCParams mpc_params = {h_, Ts_, k_hor_, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, v_min, v_max, a_min, a_max}};
        FoVCBFParams fov_cbf_params = {fov_beta, fov_Ds, fov_Rs};

        // subs
        state_sub_ = nh_.subscribe<nav_msgs::Odometry>("mavros/local_position/odom", 1, std::bind(&ControlNode::state_update_callback, this, std::placeholders::_1));
        target_subs_.resize(NUM_TARGETS);
        for (size_t i = 0; i < NUM_TARGETS; ++i) {
            size_t target_id = neighbor_ids[i];
            target_subs_[i] = nh_.subscribe<geometry_msgs::PoseWithCovarianceStamped>("target_"+std::to_string(target_id)+"/estimate", 1, std::bind(&ControlNode::target_update_callback, this, std::placeholders::_1, i));
//            target_subs_[i] = nh_.subscribe<nav_msgs::Odometry>("/uav"+std::to_string(target_id)+"/mavros/local_position/odom", 1, std::bind(&ControlNode::target_odom_callback, this, std::placeholders::_1, i));
        }
        goal_sub_ = nh_.subscribe<geometry_msgs::Pose>("goal", 1, std::bind(&ControlNode::goal_update_callback, this, std::placeholders::_1));

        // control pub
        control_pub_ = nh_.advertise<mavros_msgs::PositionTarget>("mavros/setpoint_raw/local", 10);

        // planned path visualization
        path_pub_ = nh_.advertise<nav_msgs::Path>("planned_path", 10);

        optimizer_ = nh_.createTimer(ros::Duration(h_), std::bind(&ControlNode::optimization_callback, this));
        timer_ = nh_.createTimer(ros::Duration(Ts_), std::bind(&ControlNode::timer_callback, this));
        take_off_module_ = nh_.createTimer(ros::Duration(Ts_), std::bind(&ControlNode::takeoff_callback, this));
        std::cout << "finished all sub/pub" << "\n";

        std::cout << "Ts: " << Ts_ << ", h: " << h_ << std::endl;

        // json for record
        std::string JSON_FILENAME = "../../../tools/CBFXYYawStates.json";
        json states;
        states["dt"] = h_;
        states["Ts"] = Ts_;
        // init model
        std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h_);
        std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts_);
        StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h_/Ts_));
        StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h_/Ts_));
        // init cbf
        std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, v_min, v_max);
        // init bezier mpc-cbf
        uint64_t bezier_continuity_upto_degree = 7;
        int impc_iter = 2;
        int cbf_horizon = 2;
        double slack_cost = 1000;
        double slack_decay_rate = 0.2;
        bool slack_mode = false;
        BezierMPCCBFParams bezier_mpc_cbf_params = {piecewise_bezier_params, mpc_params, fov_cbf_params};
        IMPCParams impc_params = {cbf_horizon, impc_iter, slack_cost, slack_decay_rate, slack_mode};
        IMPCCBFParams impc_cbf_params = {bezier_mpc_cbf_params, impc_params};
        bezier_impc_cbf_ptr_ = std::make_unique<BezierIMPCCBF>(impc_cbf_params, pred_model_ptr, fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr, NUM_TARGETS);
        std::cout << "successfully init the control core..." << "\n";

        last_request_ = ros::Time::now();
//        std::cout << "in init, last request time: " << last_request_.toSec() << "\n";
//        std::cout << "in init, now time: " << ros::Time::now().toSec() << "\n";
    }

    ~ControlNode()
    {
        std::cout << "ControlNode destroyer called\n";
    }

    void stop();
    void takeoff_callback();
    void timer_callback();
    void optimization_callback();
    void target_update_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& msg, size_t target_index);
    void target_odom_callback(const nav_msgs::Odometry::ConstPtr& msg, size_t target_index);
    void state_update_callback(const nav_msgs::Odometry::ConstPtr& msg);
    void goal_update_callback(const geometry_msgs::Pose::ConstPtr& msg);
    VectorDIM convertToClosestYaw();

private:
    int ROBOT_ID = 0;
    int NUM_TARGETS = 1;
    double rate = 10.0;
    std::string CONFIG_FILENAME;

    // mpc settings
    double h_;
    double Ts_;
    int k_hor_;


    geometry_msgs::PoseStamped takeoff_pose_;
    geometry_msgs::PoseStamped land_pose_;
    State state_;
    std::vector<VectorDIM> target_states_;
    std::vector<Matrix> target_covs_;
    VectorDIM goal_;
    std::shared_ptr<SingleParameterPiecewiseCurve> curve_;
    double eval_t_;
    double z_ = 1;
    std::unique_ptr<BezierIMPCCBF> bezier_impc_cbf_ptr_;
    bool taken_off_ = false;
    bool land_ = false;
    bool optim_done_ = false;
    bool last_request_time_init = false;

    // ROS
    ros::NodeHandle nh_;
    ros::NodeHandle nh_priv_;
    ros::Subscriber state_sub_;
    ros::Subscriber mavros_state_sub_;
    ros::Subscriber pose_sub_;
    std::vector<ros::Subscriber> target_subs_;
    ros::Subscriber goal_sub_;
    ros::Publisher control_pub_;
    ros::Publisher local_pos_pub_;
    ros::Publisher path_pub_;
    ros::Timer optimizer_;
    ros::Timer timer_;
    ros::Timer take_off_module_;
    ros::ServiceClient arming_client_;
    ros::ServiceClient set_mode_client_;
    ros::Time last_request_;
    ros::Time mission_start_time_;
    ros::Time land_start_time_;
    ros::Time last_traj_optim_t_;
    mavros_msgs::SetMode offb_set_mode_;

};

void ControlNode::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
}

void ControlNode::target_update_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& msg, size_t target_index) {
    target_states_[target_index] << msg->pose.pose.position.x, msg->pose.pose.position.y, 0;
    Matrix target_cov(3,3);
    target_cov << msg->pose.covariance[0], msg->pose.covariance[1], msg->pose.covariance[5],
            msg->pose.covariance[6], msg->pose.covariance[7], msg->pose.covariance[11],
            msg->pose.covariance[30], msg->pose.covariance[31], msg->pose.covariance[35];
    target_covs_[target_index] = target_cov;
}

void ControlNode::target_odom_callback(const nav_msgs::Odometry::ConstPtr& msg, size_t target_index) {
    target_states_[target_index] << msg->pose.pose.position.x, msg->pose.pose.position.y, 0;
    Matrix target_cov(3,3);
    target_cov << msg->pose.covariance[0], msg->pose.covariance[1], msg->pose.covariance[5],
                  msg->pose.covariance[6], msg->pose.covariance[7], msg->pose.covariance[11],
                  msg->pose.covariance[30], msg->pose.covariance[31], msg->pose.covariance[35];
    target_covs_[target_index] = target_cov;
}

void ControlNode::state_update_callback(const nav_msgs::Odometry::ConstPtr& msg) {
    // convert the orientation to yaw
    tf2::Quaternion q(msg->pose.pose.orientation.x, msg->pose.pose.orientation.y, msg->pose.pose.orientation.z, msg->pose.pose.orientation.w);
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);
    // update the position
    state_.pos_ << msg->pose.pose.position.x, msg->pose.pose.position.y, yaw;
    // std::cout << "Pos: " << state_.pos_.transpose() << std::endl; TODO Vicon has no access to vel directly
//    state_.vel_ << msg->twist.twist.linear.x, msg->twist.twist.linear.y, msg->twist.twist.angular.z;
    // std::cout << "Vel: " << state_.vel_.transpose() << std::endl;
}

void ControlNode::goal_update_callback(const geometry_msgs::Pose::ConstPtr& msg) {
    // update the goal
    tf2::Quaternion q(msg->orientation.x, msg->orientation.y, msg->orientation.z, msg->orientation.w);
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);
    goal_ << msg->position.x, msg->position.y, yaw;
    goal_ = convertToClosestYaw();
}

void ControlNode::optimization_callback() {
    if (taken_off_) {
        bezier_impc_cbf_ptr_->resetProblem();

        Vector ref_positions(DIM * k_hor_);
        // static target reference
        ref_positions = goal_.replicate(k_hor_, 1);
//        std::cout << "goal finish\n";


//        std::cout << "goal: " << goal_.transpose() << "\n";
        // std::cout << "state pos: " << state_.pos_.transpose() << "\n";
        // std::cout << "state vel: " << state_.vel_.transpose() << "\n";
//         std::cout << "target size: " << target_states_.size() << "\n";
//         std::cout << "ROBOT_ID:" << ROBOT_ID << "target 0: " << target_states_[0].transpose() << "\n";
//         std::cout << "ROBOT_ID:" << ROBOT_ID << "target 0 cov: " << target_covs_[0] << "\n";
        // std::cout << "target 1: " << target_states_[1].transpose() << "\n";
//        std::cout << "target: " << target_states_.size() << "\n";

        std::vector <SingleParameterPiecewiseCurve> trajs;
        bool success = bezier_impc_cbf_ptr_->optimize(trajs, state_, target_states_, target_covs_, ref_positions);
//        std::cout << "ROBOT_ID:" << ROBOT_ID << "optimize finish\n";
        if (!success) {
            std::cout << "QP not success" << "\n";
        }

        // update the plan only if the traj_optim success
        if (success) {
            curve_ = std::make_shared<SingleParameterPiecewiseCurve>(std::move(trajs.back()));
            // reset the eval_t
            eval_t_ = 0;
            last_traj_optim_t_ = ros::Time::now();
            optim_done_ = true;
        }

        // Publish planned path
        if (optim_done_) {
            nav_msgs::Path path_msg;
            path_msg.header.frame_id = "map";
            path_msg.header.stamp = ros::Time::now();
            for (double t_i = 0; t_i < h_*(k_hor_-1); t_i=t_i+2*h_)
            {
                geometry_msgs::PoseStamped pose;
                VectorDIM eval = curve_->eval(t_i, 0);
                pose.pose.position.x = eval(0);
                pose.pose.position.y = eval(1);
                pose.pose.position.z = z_;
                tf2::Quaternion q;
                q.setRPY(0, 0, eval(2));
                pose.pose.orientation = tf2::toMsg(q);
                path_msg.poses.push_back(pose);
            }
            path_pub_.publish(path_msg);
        }
    }
}

void ControlNode::takeoff_callback() {
    if (!taken_off_) {
        // arming and offboard
        offb_set_mode_.request.custom_mode = "OFFBOARD";
        mavros_msgs::CommandBool arm_cmd;
        arm_cmd.request.value = true;

        if (!last_request_time_init) {
//            std::cout << "last request time: " << last_request_.toSec() << "\n";
            last_request_time_init = true;
            last_request_ = ros::Time::now();
//            std::cout << "last request time: " << last_request_.toSec() << "\n";
        }

        if (current_state.mode != "OFFBOARD" && ros::Time::now() - last_request_ > ros::Duration(2.0)) {
            if (set_mode_client_.call(offb_set_mode_) && offb_set_mode_.response.mode_sent) {
                ROS_INFO("Offboard enabled.");
            }
            last_request_ = ros::Time::now();
        } else {
            if (!current_state.armed && (ros::Time::now() - last_request_ > ros::Duration(2.0))) {
                if (arming_client_.call(arm_cmd) && arm_cmd.response.success) {
                    ROS_INFO("Vehicle armed. Start to take off");
                }
                last_request_ = ros::Time::now();
            }
        }
        // take off
        if (ros::Time::now() - last_request_ < ros::Duration(takeoff_time)) {
            local_pos_pub_.publish(takeoff_pose_);
        } else {
            std::cout << "Finish take off, execute mission..." << "\n";
            taken_off_ = true;
            mission_start_time_ = ros::Time::now();
        }
    }
}

void ControlNode::timer_callback() {
    if (taken_off_ && !optim_done_ && !land_) {
        local_pos_pub_.publish(takeoff_pose_);
    }
    else if (taken_off_ && optim_done_ && !land_) {
        // std::cout << "Inside control loop\n";
        eval_t_ = (ros::Time::now() - last_traj_optim_t_).toSec();
        // this could happen when optimization fails
        if (eval_t_ > curve_->max_parameter()) {
            eval_t_ = curve_->max_parameter();
        }
//        std::cout << "eval time: " << eval_t_ << "\n";
        std::vector<VectorDIM> evals;
        for (size_t d = 0; d <= 2; ++d) {
            // std::cout << "eval time" << eval_t_ << "\n";
            VectorDIM eval = curve_->eval(eval_t_, d);
            // std::cout << "eval : " << eval.transpose() << "\n";
            evals.push_back(eval);
        }

        // publish the control msg
        mavros_msgs::PositionTarget msg;
        msg.header.stamp = ros::Time::now();
        msg.coordinate_frame = 1;
//        std::cout << "init state: " << init_state_(0) << ", " << init_state_(1) << ", " << init_state_(2) << "\n";
        msg.position.x = evals[0](0);
        msg.position.y = evals[0](1);
        msg.position.z = z_;
        msg.velocity.x = evals[1](0);
        msg.velocity.y = evals[1](1);
        msg.velocity.z = 0;
        msg.acceleration_or_force.x = evals[2](0);
        msg.acceleration_or_force.y = evals[2](1);
        msg.acceleration_or_force.z = 0;
        msg.yaw = evals[0](2);
        msg.yaw_rate = evals[1](2);
        control_pub_.publish(msg);
        // update the higher order states
        state_.vel_ << evals[1];
        // std::cout << "Desired vel: " << evals[1].transpose() << std::endl;
        // std::cout << "Desired acc: " << evals[2].transpose() << std::endl;
        if (ros::Time::now() - mission_start_time_ > ros::Duration(mission_time)) {
            land_pose_.header.frame_id = "map";
            land_pose_.pose.position.x = state_.pos_(0);
            land_pose_.pose.position.y = state_.pos_(1);
            land_pose_.pose.position.z = 0.2;
            double yaw = state_.pos_(2);
            tf2::Quaternion q;
            q.setEuler(0, 0, yaw);
            land_pose_.pose.orientation = tf2::toMsg(q);
            land_start_time_ = ros::Time::now();
            land_ = true;
            std::cout << "Mission finished, start landing...\n";
        }
    } else if (land_) {
         offb_set_mode_.request.custom_mode = "AUTO.LAND";
         // go to near ground, ready to land
         if (ros::Time::now() - land_start_time_ < ros::Duration(land_time)) {
             local_pos_pub_.publish(land_pose_);
         }
         // auto-land
         if (ros::Time::now() - land_start_time_ >= ros::Duration(land_time) && current_state.mode != "AUTO.LAND" &&
             ros::Time::now() - last_request_ > ros::Duration(2.0)) {
             if( set_mode_client_.call(offb_set_mode_) &&
             offb_set_mode_.response.mode_sent){
                 ROS_INFO("Auto Land Triggered...");
             }
             last_request_ = ros::Time::now();
         }
    }

}

ControlNode::VectorDIM ControlNode::convertToClosestYaw() {
//    std::cout << "start get the current_yaw\n";
    double current_yaw = state_.pos_(2);
    // generate the candidate desire yaws
    ControlNode::Vector candidate_yaws(3);

//    std::cout << "start build candidate yaw\n";
//    std::cout << goal_(2) << "\n";
//    std::cout << M_PI << "\n";
    candidate_yaws << goal_(2), goal_(2) + 2 * M_PI, goal_(2) - 2 * M_PI;

//    std::cout << "compute the offset\n";
    ControlNode::Vector candidate_yaws_offset(3);
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

    ControlNode::VectorDIM converted_goal;
    converted_goal << goal_(0), goal_(1), candidate_yaws(argmin_index);
//    std::cout << "converted_goal: " << converted_goal.transpose() << "\n";
    return converted_goal;

}

/*******************************************************************************
* Main function
*******************************************************************************/
//alternatively to a global variable to have access to the method you can make STATIC the class method interested,
//but some class function may not be accessed: "this->" method cannot be used

void nodeobj_wrapper_function(int){
    ROS_WARN("signal handler function CALLED");
    node_shutdown_request = 1;
}

int main(int argc, char* argv[])
{
    ros::init(argc, argv, "control_node", ros::init_options::NoSigintHandler);
    signal(SIGINT, nodeobj_wrapper_function);

    //Controller node_controller;
    auto node_controller = std::make_shared<ControlNode>();

    while (!node_shutdown_request){
        ros::spinOnce();
    }
    node_controller->stop();

    //ros::spin();
    //do pre-shutdown tasks
    if (ros::ok())
    {
        ROS_WARN("ROS HAS NOT BEEN PROPERLY SHUTDOWN, it is being shutdown again.");
        ros::shutdown();
    }

    return 0;
}