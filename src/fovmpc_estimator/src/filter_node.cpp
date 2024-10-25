#include <vector>
#include <math.h>
#include <chrono>
#include <signal.h>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <ros/ros.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Odometry.h>
#include "particle_filter/detail/particle_filter.h"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/utils.h>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>



#define M_PI   3.14159265358979323846  /*pi*/

using namespace std::chrono_literals;
using std::placeholders::_1;

sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed



class Node
{
    public:
        Node() : nh_priv_("~"), filter()
        {
            // ---------------- ROS params -----------------
            this->nh_priv_.getParam("ROBOT_ID", ROBOT_ID);
            this->nh_priv_.getParam("TARGET_ID", TARGET_ID);
            this->nh_priv_.getParam("rate", rate);
            this->nh_priv_.getParam("STATE_SIZE", state_size);
            this->nh_priv_.getParam("PARTICLES_NUM", PARTICLES_NUM);
            this->nh_priv_.getParam("ROBOT_FOV", FOV_DEG);
            this->nh_priv_.getParam("ROBOT_RANGE", ROBOT_RANGE);

            FOV_RAD = FOV_DEG * M_PI / 180.0;



            // ---------- Subs and Pubs -------------------
            odom_sub_ = nh_.subscribe<nav_msgs::Odometry>("/supervisor/uav"+std::to_string(ROBOT_ID)+"/odom", 1, std::bind(&Node::odom_callback, this, std::placeholders::_1));
            target_sub_ = nh_.subscribe<geometry_msgs::Pose>("/robot"+std::to_string(ROBOT_ID)+"/tag_"+std::to_string(TARGET_ID)+"/pose", 1, std::bind(&Node::target_callback, this, std::placeholders::_1));
            est_pub_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("/target_"+std::to_string(TARGET_ID)+"/estimate", 10);
            timer_ = nh_.createTimer(ros::Duration(1/rate), std::bind(&Node::timer_callback, this));

            // Init params
            Eigen::VectorXd state = Eigen::VectorXd::Ones(state_size);
            Eigen::MatrixXd initCov = 1.0*Eigen::MatrixXd::Identity(state_size, state_size);
            Eigen::MatrixXd processCov = 0.1*Eigen::MatrixXd::Identity(state_size, state_size);
            Eigen::MatrixXd measCov = 0.05*Eigen::MatrixXd::Identity(state_size, state_size);

            filter.init(PARTICLES_NUM, state, initCov, processCov, measCov);

            obs_rel.resize(state_size);
            estimate.resize(state_size);
            covariance.resize(state_size, state_size);
            obs_rel.setZero();
            estimate.setZero();
            covariance.setZero();

            samples.resize(3, PARTICLES_NUM);
            samples.setZero();

            // Self-localization
            p_i.resize(7);
            p_i.setZero();
            p_i(6) = 1.0;           // make quaternion valid

            // Relative localization
            meas_rel.resize(state_size);


        } 

        ~Node()
        {
            std::cout << "Node destroyer called\n";
        }

    void stop();
    void timer_callback();
    void target_callback(const geometry_msgs::Pose::ConstPtr& msg);
    void odom_callback(const nav_msgs::Odometry::ConstPtr& msg);
    bool insideFOV(Eigen::VectorXd robot, Eigen::VectorXd target, double fov, double range);
    // std::vector<bool> insideFOV(Eigen::VectorXd robot, Eigen::MatrixXd targets, double fov, double range);

    private:
        int ROBOT_ID = 0;
        int TARGET_ID = 1;
        double FOV_DEG = 120.0;
        double FOV_RAD;
        double ROBOT_RANGE = 5.0;
        double rate = 10.0;
        double state_size = 3;
        int PARTICLES_NUM = 100;
        

        pf::ParticleFilter filter;
        Eigen::VectorXd obs_rel;
        Eigen::VectorXd estimate;
        Eigen::MatrixXd covariance;
        Eigen::MatrixXd samples;
        Eigen::VectorXd p_i;
        Eigen::VectorXd meas_rel;

        // ROS
        ros::NodeHandle nh_;
        ros::NodeHandle nh_priv_;
        ros::Timer timer_;
        ros::Publisher est_pub_;
        ros::Subscriber target_sub_;
        ros::Subscriber odom_sub_;
        geometry_msgs::PoseWithCovarianceStamped est_msg;
        geometry_msgs::TransformStamped transformStamped;
        tf2_ros::TransformBroadcaster br;

};

void Node::target_callback(const geometry_msgs::Pose::ConstPtr& msg)
{
    meas_rel(0) = msg->position.x;
    meas_rel(1) = msg->position.y;
    meas_rel(2) = msg->position.z;
    std::cout << "Rel measurement: " << meas_rel.transpose() << std::endl;
    tf2::Quaternion q(p_i(3), p_i(4), p_i(5), p_i(6));
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
        -sin(yaw), cos(yaw), 0.0,
        0.0, 0.0, 1.0;
    std::cout << "Robot position : " << p_i.head(3).transpose() << std::endl;
    Eigen::VectorXd meas_glob = p_i.head(3) + R.transpose() * meas_rel;
    std::cout << "Global position: " << meas_glob.transpose() << std::endl;

    if (insideFOV(p_i, meas_glob, FOV_RAD, ROBOT_RANGE))
    {
        filter.update(meas_glob);
    }
}

void Node::odom_callback(const nav_msgs::Odometry::ConstPtr& msg)
{
    p_i(0) = msg->pose.pose.position.x;
    p_i(1) = msg->pose.pose.position.y;
    p_i(2) = msg->pose.pose.position.z;
    p_i(3) = msg->pose.pose.orientation.x;
    p_i(4) = msg->pose.pose.orientation.y;
    p_i(5) = msg->pose.pose.orientation.z;
    p_i(6) = msg->pose.pose.orientation.w;
}


void Node::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
}

bool Node::insideFOV(Eigen::VectorXd robot, Eigen::VectorXd target, double fov, double range)
{
    tf2::Quaternion q(robot(3), robot(4), robot(5), robot(6));
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
        -sin(yaw), cos(yaw), 0.0,
        0.0, 0.0, 1.0;

    Eigen::VectorXd t_local = R * (target - robot.head(3));
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


/*std::vector<bool> insideFOV(Eigen::VectorXd robot, Eigen::MatrixXd targets, double fov, double range)
{
    tf2::Quaternion q(robot(3), robot(4), robot(5), robot(6));
    tf2::Matrix3x3 m(q);
    double roll, pitch, yaw;
    m.getRPY(roll, pitch, yaw);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
        -sin(yaw), cos(yaw), 0.0,
        0.0, 0.0, 1.0;

    Eigen::MatrixXd t_local = R * (targets.colwise() - robot.head(3));
    Eigen::VectorXd dist = t_local.colwise().norm();
    // Eigen::VectorXd angles = atan2(targets.row(1).array(), targets.row(0).array())
    Eigen::VectorXd angles = targets.row(1).array().unaryExpr([](double x) {return atan2(x, targets.row(0)[x.index()]);});

    Eigen::Array<bool, Eigen::Dynamic, 1> insides = (angles.array() <= 0.5*fov && dist.array() <= range);
    return insides;
}*/


void Node::timer_callback()
{
    // Run PF
    filter.predict();

    // Remove particles inside FoV
    Eigen::VectorXd weights = filter.getWeights();
    samples = filter.getParticles();
    for (int i = 0; i < PARTICLES_NUM; ++i)
    {
        if (insideFOV(p_i, samples.col(i), FOV_RAD, ROBOT_RANGE))
        {
            weights(i) /= 10.0;
        }
    }

    filter.setWeights(weights);

    filter.resample();
    filter.estimateState();

    // Update variables
    estimate = filter.getState();
    std::cout << "Estimated state: " << estimate.transpose() << std::endl;
    covariance = filter.getDistribution();

    // Publish estimate msg
    est_msg.header.stamp = ros::Time::now();
    est_msg.header.frame_id = "map";
    est_msg.pose.pose.position.x = estimate(0);
    est_msg.pose.pose.position.y = estimate(1);
    est_msg.pose.pose.position.z = estimate(2);
    est_msg.pose.covariance[0] = covariance(0,0);       // xx
    est_msg.pose.covariance[1] = covariance(0,1);       // xy
    est_msg.pose.covariance[2] = covariance(0,2);       // xz
    est_msg.pose.covariance[6] = covariance(1,0);       // yx
    est_msg.pose.covariance[7] = covariance(1,1);       // yy
    est_msg.pose.covariance[8] = covariance(1,2);       // yz
    est_msg.pose.covariance[12] = covariance(2,0);       // zx
    est_msg.pose.covariance[13] = covariance(2,1);       // zy
    est_msg.pose.covariance[14] = covariance(2,2);       // zz
    est_pub_.publish(est_msg);

    // Publish TF
    /*transformStamped.header.stamp = ros::Time::now();
    transformStamped.header.frame_id = "uav_"+std::to_string(ROBOT_ID);
    transformStamped.child_frame_id = "target_"+std::to_string(TARGET_ID);
    transformStamped.transform.translation.x = meas_rel(0);
    transformStamped.transform.translation.y = meas_rel(1);
    transformStamped.transform.translation.z = meas_rel(2);
    transformStamped.transform.rotation.x = 0.0;
    transformStamped.transform.rotation.y = 0.0;
    transformStamped.transform.rotation.z = 0.0;
    transformStamped.transform.rotation.w = 1.0;
    br.sendTransform(transformStamped);*/


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
    ros::init(argc, argv, "test_node", ros::init_options::NoSigintHandler);
    signal(SIGINT, nodeobj_wrapper_function);

    //Controller node_controller;
    auto node_controller = std::make_shared<Node>();

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