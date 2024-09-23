import json
import numpy
import numpy as np
import PyQt5
from PyQt5 import *
from PyQt5 import QtCore, QtSvg, QtWidgets, QtGui
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import animation, rc, rcParams

def load_states(json_filename):
    f = open(json_filename)
    data = json.load(f)
    return data

def fov_xy(pos, radian, fov_beta, fov_range, num_points=100):
    pos_x = pos[0]
    pos_y = pos[1]
    angles = np.linspace(-fov_beta/2 + radian, fov_beta/2 + radian, num_points)
    # Calculate the points on the edges of the FOV
    x = fov_range * np.cos(angles) + pos_x
    y = fov_range * np.sin(angles) + pos_y
    return x,y

def animation2D_XYYaw(traj, dt, Ts, bbox, pred_curve=None, fov_beta=None, fov_range=None, Trailing=True, PredCurve=True, save_name= "./test.mp4"):
    # animation settings
    n_agent, total_frame, _ = traj.shape
    trailing_buf_size = 100
    dir_len = 0.1  # length for orientation tick

    # load traj
    states = traj  # [n, ts, dim]
    x = traj[:, :, 0]
    y = traj[:, :, 1]
    yaw = traj[:, :, 2]
    x_v = traj[:, :, 3]
    y_v = traj[:, :, 4]
    yaw_v = traj[:, :, 5]

    # set frame scale
    x_min = np.min(x)-3
    x_max = np.max(x)+3
    y_min = np.min(y)-3
    y_max = np.max(y)+3
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    # Set equal aspect ratio
    ax.set_aspect('equal')

    # plot the init frame
    p, = ax.plot(x[0], y[0], c='dodgerblue', marker='o', markersize=5, linestyle='None', alpha=0.6)
    vel_line = [Line2D([x[0][0], dir_len * x_v[0][0] + x[0][0]],
                       [y[0][0], dir_len * y_v[0][0] + y[0][0]]) for _ in range(n_agent)]
    for i in range(n_agent):
        vel_line[i], = ax.plot([x[i][0], dir_len * x_v[i][0] + x[i][0]],
                               [y[i][0], dir_len * y_v[i][0] + y[i][0]],
                               c='dodgerblue', alpha=0.6)

    fov = []
    for i in range(n_agent):
        pos = [x[i, 0], y[i, 0]]
        fov_x, fov_y = fov_xy(pos, yaw[i, 0], fov_beta, fov_range)
        fov.append(ax.fill(np.concatenate(([pos[0]], fov_x, [pos[0]])), np.concatenate(([pos[1]], fov_y, [pos[1]])), color='lightblue', alpha=0.5))

    # Add the rectangles (boxes)
    boxes = []
    for i in range(n_agent):
        rect = plt.Rectangle((x[i][0] - bbox[0], y[i][0] - bbox[1]), 2*bbox[0], 2*bbox[1],
                             linewidth=1, edgecolor='blue', facecolor='none')
        ax.add_patch(rect)
        boxes.append(rect)

    if Trailing:
        trail = [Line2D([x[0][0]], [y[0][0]]) for i in range(n_agent)]
        for i in range(n_agent):
            trail[i], = ax.plot([x[i, 0]], [y[i, 0]], c='dodgerblue', marker=None, alpha=0.2)

    if PredCurve:
        pred = [Line2D([x[0][0]], [y[0][0]]) for i in range(n_agent)]
        for i in range(n_agent):
            pred[i], = ax.plot([x[i, 0]], [y[i, 0]], c='red', marker=None, alpha=0.8)

    # loop the frame
    def animate(ts):
        if ts % 100 == 0:
            print("Animating frame {0}/{1}".format(ts, total_frame))

        p.set_data(x[:, ts], y[:, ts])
        for i in range(n_agent):
            vel_line[i].set_data([x[i][ts], dir_len * x_v[i][ts] + x[i][ts]],
                                 [y[i][ts], dir_len * y_v[i][ts] + y[i][ts]])

        for i in range(n_agent):
            boxes[i].set_xy((x[i][ts] - bbox[0], y[i][ts] - bbox[1]))

        for i in range(n_agent):
            pos = [x[i, ts], y[i, ts]]
            fov_x, fov_y = fov_xy(pos, yaw[i, ts], fov_beta, fov_range)
            fov[i][0].set_xy(np.column_stack((np.concatenate(([pos[0]], fov_x, [pos[0]])),
                                              np.concatenate(([pos[1]], fov_y, [pos[1]])))))

        trail_offset = trailing_buf_size if ts > trailing_buf_size else ts
        if Trailing:
            for i in range(n_agent):
                trail[i].set_data(x[i, ts - trail_offset:ts], y[i, ts - trail_offset:ts])

        pred_index = int(ts//int(dt/Ts))
        if PredCurve:
            for i in range(n_agent):
                pred[i].set_data(pred_curve[i, pred_index, :, 0], pred_curve[i, pred_index, :, 1])

        return p,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, frames=total_frame, interval=Ts*1e+3, blit=True)
    anim.save(save_name, writer=animation.FFMpegWriter(fps=1/Ts))
    return anim

def plot2D_XYYaw(traj, obs_time=None, save_name="./test.jpg"):

    n_agent = traj.shape[0]
    fig = plt.figure()
    # use position data to plot the 3D trajectory
    for obs_index in range(n_agent):
        if obs_time != None:
            states = traj[obs_index, :obs_time, :]
            x = states[:, 0]
            y = states[:, 1]
            yaw = states[:, 2]
            x_v = states[:, 3]
            y_v = states[:, 4]
            yaw_v = states[:, 5]
        else:
            states = traj[obs_index, :, :]
            x = states[:, 0]
            y = states[:, 1]
            yaw = states[:, 2]
            x_v = states[:, 3]
            y_v = states[:, 4]
            yaw_v = states[:, 5]
        plt.plot(x, y, label='motion' + str(obs_index))
    plt.legend()
    plt.title('2D Trajectory')
    plt.savefig(save_name)
    plt.show()


if __name__ == "__main__":
    num_robots = len(load_states("CBFXYYawStates.json")["robots"])
    traj = np.array([load_states("CBFXYYawStates.json")["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
    dt = load_states("CBFXYYawStates.json")["dt"]
    Ts = load_states("CBFXYYawStates.json")["Ts"]
    bbox = load_states("../config/config.json")["robot_params"]["collision_shape"]["aligned_box"][:2]
    pred_curve = np.array([load_states("CBFXYYawStates.json")["robots"][str(_)]["pred_curve"] for _ in range(num_robots)])
    # plot2D_XYYaw(traj)
    FoV_beta = load_states("../config/config.json")["fov_cbf_params"]["beta"] * np.pi/180
    FoV_range = load_states("../config/config.json")["fov_cbf_params"]["Rs"]
    animation2D_XYYaw(traj, dt, Ts, bbox, pred_curve, FoV_beta, FoV_range)