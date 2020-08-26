#include <memory>
#include <iostream>
#include <vector>
#include <chrono>
#include <utility>
#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/pose.hpp>
#include <Trackgraph.h>
#include <Obstacle.h>
#include <ParseGraph.h>
#include <Path_planning_utils.h>
#include <cubic_spline.h>

using std::placeholders::_1;

class PathPlanner : public rclcpp::Node
{
  public:
    PathPlanner(Trackgraph *global_graph)
    : Node("minimal_subscriber")
    {
      subscription_ = this->create_subscription<geometry_msgs::msg::Pose>(
      "pose", 10, std::bind(&PathPlanner::planner_callback, this, _1));
      global_graph_ = global_graph;
    }
    int plan_horizon = 300; // meters

  private:
    void planner_callback(const geometry_msgs::msg::Pose::SharedPtr msg) 
    {
      Obstacles obstacles; // TODO: read in obstacles from another node
      obstacles.sv.push_back(StaticObstacle{3.75, 1.46, 1.8});

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
        
      RCLCPP_INFO(this->get_logger(), "Traj C2 X-Spline Coeffs[3]: a='%f', b = '%f', c = '%f', d = '%f'",
          c2_traj.first[3].a, c2_traj.first[3].b, c2_traj.first[3].c, c2_traj.first[3].d);
      RCLCPP_INFO(this->get_logger(), "Start Pose x: '%f', y: '%f'", msg->position.x, msg->position.y);
      RCLCPP_INFO(this->get_logger(), "End Path Node x: '%f', y: '%f'",
                  mcp_xy.first.back(), mcp_xy.second.back());
    }
    rclcpp::Subscription<geometry_msgs::msg::Pose>::SharedPtr subscription_;
    Trackgraph *global_graph_;
};

int main(int argc, char * argv[])
{
  rclcpp::init(argc, argv);

  // Track lattice input files   
  string basePath = "/home/mschoder/mit_driverless/indy_pers/indyCar_obstacle_avoidance/";
  string nodeFile = basePath + "data/nodes.json";
  string edgeFile = basePath + "data/edges.json";

  // Parse global graph json data
  Trackgraph global_graph = parseGlobalGraph(nodeFile, edgeFile);
  Trackgraph *gg_ptr = &global_graph;

  // Callback performs all local graph extraction and processing
  rclcpp::spin(std::make_shared<PathPlanner>(gg_ptr));
  rclcpp::shutdown();
  return 0;
}
