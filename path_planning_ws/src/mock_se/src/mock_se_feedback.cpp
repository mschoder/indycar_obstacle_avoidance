#include <memory>
#include <chrono>
#include <functional>

#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/pose.hpp"
#include "nav_msgs/msg/path.hpp"


using std::placeholders::_1;
using namespace std::chrono_literals;

class MockSEfeedback : public rclcpp::Node {
    public:
        MockSEfeedback(): Node("mock_se_feedback") {
            sub_ = this->create_subscription<nav_msgs::msg::Path>(
                "traj", 10, std::bind(&MockSEfeedback::topic_callback, this, _1));
                 
            pub_ = this->create_publisher<geometry_msgs::msg::Pose>("pose", 10);
            timer_ = this->create_wall_timer(
                500ms, std::bind(&MockSEfeedback::timer_callback, this));
        }

    private:
        void topic_callback(const nav_msgs::msg::Path::SharedPtr msg) {
            x_loc = msg->poses[IDX_JUMP].pose.position.x;
            y_loc = msg->poses[IDX_JUMP].pose.position.y;

            RCLCPP_INFO(this->get_logger(), 
            "Sending back pose['%f'] as new position. x: '%f', y: '%f'", IDX_JUMP,
            msg->poses[IDX_JUMP].pose.position.x, msg->poses[IDX_JUMP].pose.position.y);
        }

        void timer_callback() {
            geometry_msgs::msg::Pose msg;
            msg.position.x = x_loc;
            msg.position.y = y_loc;
            RCLCPP_INFO(this->get_logger(), "Publishing new pose: x = '%f', y = '%f'", 
                msg.position.x, msg.position.y);
            pub_->publish(msg);
        }

        rclcpp::TimerBase::SharedPtr timer_;
        rclcpp::Subscription<nav_msgs::msg::Path>::SharedPtr sub_;
        rclcpp::Publisher<geometry_msgs::msg::Pose>::SharedPtr pub_;
        double x_loc = 0;
        double y_loc = 0;
        const int IDX_JUMP =1;
};

int main(int argc, char * argv[]) {
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<MockSEfeedback>());
    rclcpp::shutdown();
    return 0;
}
