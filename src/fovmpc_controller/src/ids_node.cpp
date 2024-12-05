#include <nlohmann/json.hpp>
#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/BezierIMPCCBF.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <ros/ros.h>
#include <std_msgs/Int32MultiArray.h>
#include <signal.h>
#include <Eigen/Dense>
#include <Eigen/Core>


#define M_PI   3.14159265358979323846  /*pi*/

using namespace std::chrono_literals;
using std::placeholders::_1;

sig_atomic_t volatile node_shutdown_request = 0;    //signal manually generated when ctrl+c is pressed

class IDPub
{
    public:
        using json = nlohmann::json;
        IDPub(): nh_priv_("~")
        {
            this->nh_priv_.getParam("CONFIG_FILENAME", CONFIG_FILENAME);
            std::string experiment_config_filename = CONFIG_FILENAME;
            std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
            json experiment_config_json = json::parse(experiment_config_fc);
            std::cout << "successfully load the json..." << "\n";

            ROBOTS_NUM = experiment_config_json["vision"]["tags"].size();
            std::cout << "Number of robts: " << ROBOTS_NUM << std::endl;
            ids.resize(ROBOTS_NUM, 2);
            ids.setZero();
            pubs_.resize(ROBOTS_NUM);
            for (int i = 0; i < ROBOTS_NUM; i++)
            {
                pubs_[i] = nh_.advertise<std_msgs::Int32MultiArray>("/uav"+std::to_string(i)+"/tags", 10);
                ids(i, 0) = experiment_config_json["vision"]["tags"][i][0];
                ids(i, 1) = experiment_config_json["vision"]["tags"][i][1];
                std::cout << "ids of robot " << i<< ": " << ids.row(i);
            }

            timer_ = nh_.createTimer(ros::Duration(0.25), std::bind(&IDPub::timer_callback, this));
            
        }

        ~IDPub()
        {
            std::cout << "IDPub destroyer called\n";
        }

        void stop();
        void timer_callback();


    private:
        std::string CONFIG_FILENAME = ;
        int ROBOTS_NUM = 3;
        Eigen::MatrixXd ids;

        // ROS
        ros::NodeHandle nh_;
        ros::NodeHandle nh_priv_;
        ros::Timer timer_;
        std::vector<ros::Publisher> pubs_;
};

void IDPub::stop()
{
    ros::Duration(0.1).sleep();
    ros::shutdown();
}


void IDPub::timer_callback()
{
    for (int i = 0; i < ROBOTS_NUM; i++)
    {
        std_msgs::Int32MultiArray msg;
        msg.layout.dim.resize(1);
        msg.layout.dim[0].size = 2;
        for (int j = 0; j < ids.cols(); j++)
        {
            msg.data.push_back(ids(i,j));
        }
        pubs_[i].publish(msg);
    }

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
    ros::init(argc, argv, "IDPub_node", ros::init_options::NoSigintHandler);
    signal(SIGINT, nodeobj_wrapper_function);

    //Controller node_controller;
    auto node_controller = std::make_shared<IDPub>();

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