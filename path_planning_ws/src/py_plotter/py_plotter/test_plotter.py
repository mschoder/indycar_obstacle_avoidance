
import rclpy
import sys
import os
import pathlib

from rclpy.node import Node
from std_msgs.msg import String
from geometry_msgs.msg import Pose
from nav_msgs.msg import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
sys.path.append('../../py_plotter')
from py_plotter import plot_utils

## weird hack to get to root directory... TODO-fix
track_fpath = os.path.join(pathlib.Path(__file__).parents[6], 'data', 'ims_track_parsed.csv')


class MinimalSubscriber(Node):

    def __init__(self):
        super().__init__('test_plotter')
        self.subscription = self.create_subscription(
            Pose,
            '/pose',
            self.listener_callback,
            10)
        self.subscription  # prevent unused variable warning

        self.traj_sub = self.create_subscription(Path, '/traj', self.traj_cb, 10)
        self.traj_sub

        td = plot_utils.load_parsed_trackdata(track_fpath)
        ax1 = plot_utils.plot_trackfile(td)
        self.ax = ax1
        self.gx, self.gy = [], []
        self.tx, self.ty = [], []
        self.plot, = self.ax.plot(self.gx, self.gy, 'ro')
        self.plot_traj, = self.ax.plot(self.tx, self.ty, 'go')
        self.plot_traj2, = self.ax.plot(self.tx, self.ty, 'b-')
        plt.ion()
        plt.draw()

    def listener_callback(self, msg):

        # self.get_logger().info('I heard: "%s"' % msg.position.x)

        x_pos = msg.position.x
        y_pos = msg.position.y
        self.plot.set_data(x_pos, y_pos)
        plt.draw()
        plt.pause(0.00000000001)

    def traj_cb(self, msg):
        tx, ty = [], []
        for pose in msg.poses:
            tx.append(pose.pose.position.x)
            ty.append(pose.pose.position.y)
        self.plot_traj.set_data(tx, ty)
        self.plot_traj2.set_data(tx, ty)
        plt.draw()
        plt.pause(0.00000000001)
            


def main(args=None):
    rclpy.init(args=args)

    minimal_subscriber = MinimalSubscriber()
    rclpy.spin(minimal_subscriber)

    minimal_subscriber.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()