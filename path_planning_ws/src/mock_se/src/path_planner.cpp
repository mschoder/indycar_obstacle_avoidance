#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <utility>
#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/pose.hpp>
#include <trajectory_msgs/msg/joint_trajectory.hpp>
#include <nav_msgs/msg/path.hpp>

#include <Trackgraph.h>
#include <Obstacle.h>
#include <ParseGraph.h>
#include <Path_planning_utils.h>
#include <cubic_spline.h>

using std::placeholders::_1;
rclcpp::Clock Clock;

class PathPlanner : public rclcpp::Node {
    public:
        PathPlanner(Trackgraph *global_graph): Node("minimal_subscriber") {
            subscription_ = this->create_subscription<geometry_msgs::msg::Pose>(
                "pose", 10, std::bind(&PathPlanner::planner_callback, this, _1));  

            publisher_ = this->create_publisher<nav_msgs::msg::Path>("traj", 10);
            global_graph_ = global_graph;
        }
        int plan_horizon = 300; // meters

    private:
        void planner_callback(const geometry_msgs::msg::Pose::SharedPtr msg) {
            Obstacles obstacles; // TODO: read in obstacles from another node
            obstacles.sv.push_back(StaticObstacle{1.5, -150, 3.5});
            obstacles.sv.push_back(StaticObstacle{-1.0, -220, 5.0});
            obstacles.sv.push_back(StaticObstacle{7.0, -260, 5.0});
            obstacles.sv.push_back(StaticObstacle{80, -480, 5.0});
            obstacles.sv.push_back(StaticObstacle{150, -520, 5.0});
            obstacles.sv.push_back(StaticObstacle{210, -550, 5.0});

            double cur_x = msg->position.x;
            double cur_y = msg->position.y;
            pair<int, double> cur_frenet = global_graph_->xy2frenet(cur_x, cur_y);
            Trackgraph local_graph = global_graph_->extractLocalGraph(
                cur_frenet.first, cur_frenet.second, plan_horizon);

            pp_utils::collision_checker(local_graph, obstacles);
            vector<pair<int, double>> mcp = local_graph.min_cost_path_search();
            pair<vector<double>, vector<double>> mcp_xy = mcp_to_xy(local_graph, mcp);

            pair<vector<SplineSet>, vector<SplineSet>> c2_traj = pp_utils::trajectory_gen(
                local_graph, mcp);
        
            /* Debug logging */
            RCLCPP_INFO(this->get_logger(), "Traj C2 X-Spline Coeffs[3]: a='%f', b = '%f', c = '%f', d = '%f'",
                c2_traj.first[3].a, c2_traj.first[3].b, c2_traj.first[3].c, c2_traj.first[3].d);
            RCLCPP_INFO(this->get_logger(), "Start Pose x: '%f', y: '%f'", msg->position.x, msg->position.y);
            RCLCPP_INFO(this->get_logger(), "End Path Node x: '%f', y: '%f'",
                        mcp_xy.first.back(), mcp_xy.second.back());

            // Create trajectory message & publish (not using C2 spline info though)
            nav_msgs::msg::Path traj;
            geometry_msgs::msg::PoseStamped pose;
            for (int i = 0; i < mcp_xy.first.size(); ++i) {
                pose.header.frame_id = to_string(i);
                pose.pose.position.x = mcp_xy.first[i];
                pose.pose.position.y = mcp_xy.second[i];
                traj.poses.push_back(pose);
            } 
            // traj.header.stamp = 123123;
            publisher_->publish(traj);
    

            //tmp eval - rough, using original points
            // vector<double> xspl, yspl;
            // for (int i = 0; i < c2_traj.first.size(); ++i) {
            //     xspl.push_back(c2_traj.first[i].x)
            // }

            }

        rclcpp::Subscription<geometry_msgs::msg::Pose>::SharedPtr subscription_;
        rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr publisher_;
        Trackgraph *global_graph_;
};

int main(int argc, char * argv[]) {
    

    // Track lattice input files   
    string basePath = "/home/mschoder/mit_driverless/indy_pers/indyCar_obstacle_avoidance/path_planning_ws/";
    string nodeFile = basePath + "data/nodes.json";
    string edgeFile = basePath + "data/edges.json";

    rclcpp::init(argc, argv);

    // Parse global graph json data
    Trackgraph global_graph = parseGlobalGraph(nodeFile, edgeFile);
    Trackgraph *gg_ptr = &global_graph;

    // Callback performs all local graph extraction and processing
    rclcpp::spin(std::make_shared<PathPlanner>(gg_ptr));
    rclcpp::shutdown();
    return 0;
}
