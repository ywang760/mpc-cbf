//
// Created by lishuo on 10/21/24.
//

#include <math/Types.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

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
#include <signal.h>

using std::placeholders::_1;
sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed

constexpr unsigned int DIM = 3U;
class GoalNode
{
public:
    using VectorDIM = math::VectorDIM<double, DIM>;
    using json = nlohmann::json;

    GoalNode() : nh_priv_("~")
    {
        // ---------------- ROS params -----------------
        this->nh_priv_.getParam("ROBOT_ID", ROBOT_ID);
        this->nh_priv_.getParam("rate", rate);
        this->nh_priv_.getParam("CONFIG_FILENAME", CONFIG_FILENAME);


        // ---------- Subs and Pubs -------------------
        // control pub
        goal_pub_ = nh_.advertise<geometry_msgs::Pose>("/uav" + std::to_string(ROBOT_ID) + "/goal", 10);
        timer_ = nh_.createTimer(ros::Duration(1/rate), std::bind(&GoalNode::timer_callback, this));

        // Read the goal
        // load experiment config
        std::string experiment_config_filename = CONFIG_FILENAME;
        std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
        json experiment_config_json = json::parse(experiment_config_fc);

        json sf_json = experiment_config_json["tasks"]["sf"];
        goal_ << sf_json[ROBOT_ID][0], sf_json[ROBOT_ID][1], sf_json[ROBOT_ID][2];
    }

    ~GoalNode()
    {
        std::cout << "GoalNode destroyer called\n";
    }

    void stop();
    void timer_callback();

private:
    int ROBOT_ID = 0;
    double rate = 10.0;
    std::string CONFIG_FILENAME;

    VectorDIM goal_;
    double z_ = 1;

    // ROS
    ros::NodeHandle nh_;
    ros::NodeHandle nh_priv_;
    ros::Publisher goal_pub_;
    ros::Timer timer_;
};

void GoalNode::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
}

void GoalNode::timer_callback() {
    tf2::Quaternion q;
    q.setRPY(0, 0, goal_(2));

    geometry_msgs::Pose msg;
    msg.position.x = goal_(0);
    msg.position.y = goal_(1);
    msg.position.z = z_;

    msg.orientation.x = q.getX();
    msg.orientation.y = q.getY();
    msg.orientation.z = q.getZ();
    msg.orientation.w = q.getW();
    goal_pub_.publish(msg);
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
    ros::init(argc, argv, "goal_node", ros::init_options::NoSigintHandler);
    signal(SIGINT, nodeobj_wrapper_function);

    //Controller node_controller;
    auto node_controller = std::make_shared<GoalNode>();

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