#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <fstream>

#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/pose.hpp"
// #include "tf/tf.h"
// #include "nav_msgs/msg/odometry.hpp"
#include "csv_reader.h"

using namespace std::chrono_literals;

class MockSEPublisher : public rclcpp::Node
{
  public:
    MockSEPublisher(std::string filename)
    : Node("minimal_publisher"), idx_(0)
    {
      filename_ = filename;
      read_trackdata();
      publisher_ = this->create_publisher<geometry_msgs::msg::Pose>("pose", 10);
      timer_ = this->create_wall_timer(
        500ms, std::bind(&MockSEPublisher::timer_callback, this));
    }

  private:
    
    void read_trackdata() 
    {
      CSVReader reader(filename_);
      std::vector<std::vector<double> > dataList = reader.getData();
      for (std::vector<double> vec : dataList) { // rows
        s_.push_back(vec[0]);
        x_.push_back(vec[1]);
        y_.push_back(vec[2]);
        psi_.push_back(vec[3]);
      }
    }

    void timer_callback()
    {
      auto message = geometry_msgs::msg::Pose();
      message.position.x = x_[idx_];
      message.position.y = y_[idx_];
      message.orientation.z = psi_[idx_];
      RCLCPP_INFO(this->get_logger(), "Publishing: x = '%f', y = '%f'", 
            message.position.x, message.position.y);
      publisher_->publish(message);
      idx_ += 20;
      if (idx_ >= x_.size()-1) {
          idx_ = idx_ - x_.size() - 1;
      } 
    }

    rclcpp::TimerBase::SharedPtr timer_;
    rclcpp::Publisher<geometry_msgs::msg::Pose>::SharedPtr publisher_;
    size_t idx_;
    std::vector<double> s_, x_, y_, psi_;
    std::string filename_;
  };

  int main(int argc, char * argv[])
  {
    std::string trackfile = "/home/mschoder/data/ims_track_parsed.csv";
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<MockSEPublisher>(trackfile));
    rclcpp::shutdown();
    return 0;
  }