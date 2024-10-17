#include <vector>
#include <math.h>
#include <chrono>
#include <signal.h>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <ros/ros.h>
#include "particle_filter/detail/particle_filter.h"



#define M_PI   3.14159265358979323846  /*pi*/

using namespace std::chrono_literals;
using std::placeholders::_1;

sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed



class Node
{
    public:
        Node() : nh_priv_("~")
        {
            ROS_INFO("CIAO");    
            int samples_num = 100;
            Eigen::VectorXd state = Eigen::Vector3d::Ones();
            std::cout << "Init state: " << state << std::endl;
            Eigen::MatrixXd initCov = 1.0*Eigen::Matrix3d::Identity();
            Eigen::MatrixXd processCov = 1.0*Eigen::Matrix3d::Identity();
            Eigen::MatrixXd measCov = 0.5*Eigen::Matrix3d::Identity();

            pf::ParticleFilter filter(samples_num, state, initCov, processCov, measCov);


        } 

        ~Node()
        {
            std::cout << "Node destroyer called\n";
        }

    void stop();


    private:
        int ROBOTS_NUM = 3;

        // ROS
        ros::NodeHandle nh_;
        ros::NodeHandle nh_priv_;
        ros::Timer timer_;

};

void Node::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
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