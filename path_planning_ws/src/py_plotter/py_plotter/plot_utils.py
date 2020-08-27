import numpy as np
import matplotlib.pyplot as plt


def load_parsed_trackdata(track_file):
    # data in meters; reads csv file
    track_data = np.genfromtxt(track_file, delimiter=",", names=True)
    return track_data

def plot_trackfile(track_data):
    ax = plt.gca()

    # partial zoom in of corner
    ax.set_xlim((-50, 400))
    ax.set_ylim((-600, 0))

    # ax.plot(track_data["x_mid"], track_data["y_mid"], 'b--')
    ax.plot(track_data["x_left"], track_data["y_left"], 'k-')
    ax.plot(track_data["x_right"], track_data["y_right"], 'k-')

    # plot some obstacles manually
    obs_list = [(1.5, -150, 3.5), (-1.0, -220, 5.0), (7.0, -260, 5.0), 
        (80, -480, 5.0), (150, -520, 5.0), (210, -550, 5.0)]
    plot_obstacles(ax, obs_list)

    ax.set_aspect('equal')
    return ax


def plot_obstacles(ax, obs_list):
    for obs in obs_list:
        obs_shape = plt.Circle((obs[0], obs[1]), obs[2], color = 'brown', zorder=10)
        ax.add_artist(obs_shape)


# fpath = "../../../data/ims_track_parsed.csv"
# td = load_parsed_trackdata(fpath)
# ax1 = plot_trackfile(td)
# plt.show()