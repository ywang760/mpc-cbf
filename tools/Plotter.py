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

def animation2D_XYYaw(traj, dt, pred_curve, Trailing=True, PredCurve=True, save_name= "./test.mp4"):
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

    # plot the init frame
    p, = ax.plot(x[0], y[0], c='dodgerblue', marker='o', markersize=5, linestyle='None', alpha=0.6)
    vel_line = [Line2D([x[0][0], dir_len * x_v[0][0] + x[0][0]],
                       [y[0][0], dir_len * y_v[0][0] + y[0][0]]) for _ in range(n_agent)]
    for i in range(n_agent):
        vel_line[i], = ax.plot([x[i][0], dir_len * x_v[i][0] + x[i][0]],
                               [y[i][0], dir_len * y_v[i][0] + y[i][0]],
                               c='dodgerblue', alpha=0.6)

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

        trail_offset = trailing_buf_size if ts > trailing_buf_size else ts
        if Trailing:
            for i in range(n_agent):
                trail[i].set_data(x[i, ts - trail_offset:ts], y[i, ts - trail_offset:ts])

        pred_index = int(ts//10)
        if PredCurve:
            for i in range(n_agent):
                pred[i].set_data(pred_curve[i, pred_index, :, 0], pred_curve[i, pred_index, :, 1])

        return p,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, frames=total_frame, interval=dt*1e+3, blit=True)
    anim.save(save_name, writer=animation.FFMpegWriter(fps=1/dt))
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
    num_robots = len(load_states("XYYawStates.json")["robots"])
    traj = np.array([load_states("XYYawStates.json")["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
    dt = load_states("XYYawStates.json")["dt"]
    pred_curve = np.array([load_states("XYYawStates.json")["robots"][str(_)]["pred_curve"] for _ in range(num_robots)])
    # plot2D_XYYaw(traj)
    animation2D_XYYaw(traj, dt, pred_curve)
    print()