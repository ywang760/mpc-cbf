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
#include <nav_msgs/Odometry.h>
#include <tf2/utils.h>
#include <tf2/LinearMath/Quaternion.h>

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <math.h>
#include <chrono>
#include <signal.h>

#define M_PI   3.14159265358979323846  /*pi*/

using namespace std::chrono_literals;
using std::placeholders::_1;
sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed

constexpr unsigned int DIM = 3U;
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
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    ControlNode() : nh_priv_("~")
    {
        // ---------------- ROS params -----------------
        this->nh_priv_.getParam("ROBOT_ID", ROBOT_ID);
        this->nh_priv_.getParam("NUM_TARGETS", NUM_TARGETS);
        this->nh_priv_.getParam("rate", rate);


        // ---------- Subs and Pubs -------------------
        std::vector<size_t> neighbor_ids;
        for (size_t i = 0; i < NUM_TARGETS+1; ++i) {
            if (i == ROBOT_ID) {
                continue;
            }
            neighbor_ids.push_back(i);
        }
        assert(neighbor_ids.size() == NUM_TARGETS);

        // subs
        state_sub_ = nh_.subscribe<nav_msgs::Odometry>("/supervisor/uav" + std::to_string(ROBOT_ID) + "/odom", 1, std::bind(&ControlNode::state_update_callback, this, std::placeholders::_1));
        target_subs_.resize(NUM_TARGETS);
        for (size_t i = 0; i < NUM_TARGETS; ++i) {
            size_t target_id = neighbor_ids[i];
            target_subs_[i] = nh_.subscribe<geometry_msgs::PoseWithCovarianceStamped>("/target"+std::to_string(target_id)+"estimate", 1, std::bind(&ControlNode::target_update_callback, this, std::placeholders::_1, i));
        }
        goal_sub_ = nh_.subscribe<geometry_msgs::Pose>("/uav" + std::to_string(ROBOT_ID) + "/goal", 1, std::bind(&ControlNode::goal_update_callback, this, std::placeholders::_1));

        // control pub
        control_pub_ = nh_.advertise<mavros_msgs::PositionTarget>("mavros/setpoint_raw/local", 10);

        optimizer_ = nh_.createTimer(ros::Duration(h_), std::bind(&ControlNode::optimization_call_back, this));
        timer_ = nh_.createTimer(ros::Duration(Ts_), std::bind(&ControlNode::timer_callback, this));

        // Init params
        state_.vel_ = VectorDIM::Zero();
        // load experiment config
        std::string experiment_config_filename = "../../../config/config.json";
        std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
        json experiment_config_json = json::parse(experiment_config_fc);

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
        double fov_Ds = experiment_config_json["fov_cbf_params"]["Ds"];
        double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

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
        // create params
        PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
        MPCParams mpc_params = {h_, Ts_, k_hor_, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, a_min, a_max}};
        FoVCBFParams fov_cbf_params = {fov_beta, fov_Ds, fov_Rs};

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
        std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs);
        // init bezier mpc-cbf
        uint64_t bezier_continuity_upto_degree = 4;
        int impc_iter = 2;
        BezierMPCCBFParams bezier_impc_cbf_params = {piecewise_bezier_params, mpc_params, fov_cbf_params};

        bezier_impc_cbf_ptr_ = std::make_unique<BezierIMPCCBF>(bezier_impc_cbf_params, pred_model_ptr, fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr, impc_iter);
    }

    ~ControlNode()
    {
        std::cout << "ControlNode destroyer called\n";
    }

    void stop();
    void timer_callback();
    void optimization_call_back();
    void target_update_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& msg, size_t target_index);
    void state_update_callback(const geometry_msgs::PoseStamped::ConstPtr& msg);
    void goal_update_callback(const geometry_msgs::Pose::ConstPtr& msg);

private:
    int ROBOT_ID = 0;
    int NUM_TARGETS = 0;
    double rate = 10.0;

    // mpc settings
    double h_;
    double Ts_;
    int k_hor_;

    State state_;
    std::vector<VectorDIM> target_states_;
    VectorDIM goal_;
    std::unique_ptr<SingleParameterPiecewiseCurve> curve_;
    double eval_t_;
    double z_;
    std::unique_ptr<BezierIMPCCBF> bezier_impc_cbf_ptr_;

    // ROS
    ros::NodeHandle nh_;
    ros::NodeHandle nh_priv_;
    ros::Subscriber state_sub_;
    std::vector<ros::Subscriber> target_subs_;
    ros::Subscriber goal_sub_;
    ros::Publisher control_pub_;
    ros::Timer optimizer_;
    ros::Timer timer_;

};

void ControlNode::target_update_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& msg, size_t target_index) {
    target_states_[target_index] << msg->pose.pose.position.x, msg->pose.pose.position.y, 0;
}

void ControlNode::state_update_callback(const geometry_msgs::PoseStamped::ConstPtr& msg) {
    // convert the orientation to yaw
    tf2::Quaternion q(msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z, msg->pose.orientation.w);
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);
    // update the position
    state_.pos_ << msg->pose.position.x, msg->pose.position.y, yaw;
}

void ControlNode::goal_update_callback(const geometry_msgs::Pose::ConstPtr& msg) {
    // update the goal TODO update yaw
    goal_ << msg->position.x, msg->position.y, 0;
}

void ControlNode::optimization_call_back() {
    bezier_impc_cbf_ptr_->resetProblem();

    // TODO modify the init_states & add the subscribe of ref_position.
    Vector ref_positions(DIM*k_hor_);
    // static target reference
    ref_positions = goal_.replicate(k_hor_, 1);

    std::vector<SingleParameterPiecewiseCurve> trajs;
    bool success = bezier_impc_cbf_ptr_->optimize(trajs, state_, target_states_, ref_positions);
    curve_ = std::make_unique<SingleParameterPiecewiseCurve>(std::move(trajs.back()));

    // reset the eval_t
    eval_t_ = 0;
}

void ControlNode::timer_callback() {
    std::vector<VectorDIM> evals;
    for (size_t d = 0; d < 2; ++d) {
        VectorDIM eval = curve_->eval(eval_t_, d);
        evals.push_back(eval);
    }

    // publish the control msg
    mavros_msgs::PositionTarget msg;
    msg.header.stamp = ros::Time::now();
    msg.coordinate_frame = 1;
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

    eval_t_ += Ts_;
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
    ros::init(argc, argv, "supervisor_node", ros::init_options::NoSigintHandler);
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