#include <memory>
#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/pose.hpp"


using std::placeholders::_1;

class MockSESubscriber : public rclcpp::Node
{
  public:
    MockSESubscriber()
    : Node("minimal_subscriber")
    {
      subscription_ = this->create_subscription<geometry_msgs::msg::Pose>(
      "pose", 10, std::bind(&MockSESubscriber::topic_callback, this, _1));
    }

  private:
    void topic_callback(const geometry_msgs::msg::Pose::SharedPtr msg) const
    {
      RCLCPP_INFO(this->get_logger(), "Reading pose x: '%f', y: '%f'", msg->position.x, msg->position.y);
    }
    rclcpp::Subscription<geometry_msgs::msg::Pose>::SharedPtr subscription_;
};

int main(int argc, char * argv[])
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<MockSESubscriber>());
  rclcpp::shutdown();
  return 0;
}
