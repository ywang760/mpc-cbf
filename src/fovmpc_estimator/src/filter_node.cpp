#include <vector>
#include <math.h>
#include <chrono>
#include <signal.h>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <ros/ros.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include "particle_filter/detail/particle_filter.h"



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

            // ---------- Subs and Pubs -------------------
            target_sub_ = nh_.subscribe<geometry_msgs::Pose>("/robot"+std::to_string(ROBOT_ID)+"/tag"+std::to_string(TARGET_ID)+"/pose", 1, std::bind(&Node::target_callback, this, std::placeholders::_1));
            est_pub_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("/target"+std::to_string(TARGET_ID)+"estimate", 10);
            timer_ = nh_.createTimer(ros::Duration(1/rate), std::bind(&Node::timer_callback, this));

            // Init params
            int samples_num = 100;
            Eigen::VectorXd state = Eigen::Vector2d::Ones();
            Eigen::MatrixXd initCov = 1.0*Eigen::Matrix2d::Identity();
            Eigen::MatrixXd processCov = 1.0*Eigen::Matrix2d::Identity();
            Eigen::MatrixXd measCov = 0.5*Eigen::Matrix2d::Identity();

            pf::ParticleFilter filter(samples_num, state, initCov, processCov, measCov);

            obs_rel.setZero();
            estimate.setZero();
            covariance.setZero();


        } 

        ~Node()
        {
            std::cout << "Node destroyer called\n";
        }

    void stop();
    void timer_callback();
    void target_callback(const geometry_msgs::Pose::ConstPtr& msg);


    private:
        int ROBOT_ID = 0;
        int TARGET_ID = 1;
        double rate = 10.0;

        pf::ParticleFilter filter;
        Eigen::Vector2d obs_rel;
        Eigen::Vector2d estimate;
        Eigen::Matrix2d covariance;

        // ROS
        ros::NodeHandle nh_;
        ros::NodeHandle nh_priv_;
        ros::Timer timer_;
        ros::Publisher est_pub_;
        ros::Subscriber target_sub_;
        geometry_msgs::PoseWithCovarianceStamped est_msg;

};

void Node::target_callback(const geometry_msgs::Pose::ConstPtr& msg)
{
    Eigen::Vector2d meas;
    meas << msg->position.x, msg->position.y;
    filter.update(meas);
}

void Node::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
}

void Node::timer_callback()
{
    // Run PF
    filter.predict();
    filter.resample();
    filter.estimateState();

    // Update variables
    estimate = filter.getState();
    covariance = filter.getDistribution();
    std::cout << "Estimate state: " << estimate.transpose() << std::endl;
    std::cout << "Covariance: " << covariance << std::endl;

    // Publish estimate msg
    est_msg.header.stamp = ros::Time::now();
    est_msg.header.frame_id = "target"+std::to_string(TARGET_ID);
    est_msg.pose.pose.position.x = estimate(0);
    est_msg.pose.pose.position.y = estimate(1);
    est_msg.pose.covariance[0] = covariance(0,0);       // xx
    est_msg.pose.covariance[1] = covariance(0,1);       // xy
    est_msg.pose.covariance[6] = covariance(1,0);       // yx
    est_msg.pose.covariance[7] = covariance(1,1);       // yy
    est_pub_.publish(est_msg);

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