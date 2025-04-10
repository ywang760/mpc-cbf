import json
import numpy as np
# import matplotlib
# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Circle
from matplotlib import animation
import colorsys
import argparse
import os

# visualization parameters
preview_fov_range = 0.5
preview_skip_step = 5
enable_preview = True

def generate_rgb_colors(num_colors):
    output = []
    num_colors += 1 # to avoid the first color
    for index in range(1, num_colors):
        incremented_value = 1.0 * index / num_colors
        output.append(colorsys.hsv_to_rgb(incremented_value, 0.75, 0.75))
    return np.asarray(output)

def load_states(json_filename):
    f = open(json_filename)
    data = json.load(f)
    return data

def compute_ellipse(pos, cov, nstd=2):
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    return pos, width, height, theta

def plot_cov_ellipse_2d(pos, cov, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        cov : The 2x2 covariance matrix to base the ellipse on
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    if ax is None:
        ax = plt.gca()

    pos, width, height, theta = compute_ellipse(pos, cov)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip

def fov_xy(pos, radian, fov_beta, fov_range, num_points=100):
    pos_x = pos[0]
    pos_y = pos[1]
    angles = np.linspace(-fov_beta/2 + radian, fov_beta/2 + radian, num_points)
    # Calculate the points on the edges of the FOV
    x = fov_range * np.cos(angles) + pos_x
    y = fov_range * np.sin(angles) + pos_y
    return x,y

def snapshots2D_XYYaw(traj, goals, goal_radius, estimate_mean, estimate_cov, p_near, dt, Ts, bbox, frame_idx, pred_curve=None, fov_beta=None, fov_range=None, Trailing=True, Estimation=True, PredCurve=True, ShowGoals=False, colors="r", pre_save_name= "./test"):
    n_agent, total_frame, _ = traj.shape
    # load traj
    x = traj[:, :, 0]
    y = traj[:, :, 1]
    yaw = traj[:, :, 2]
    x_v = traj[:, :, 3]
    y_v = traj[:, :, 4]
    yaw_v = traj[:, :, 5]
    if Estimation:
        est_x = estimate_mean[:, :, :, 0]  # [n, n-1, ts]
        est_y = estimate_mean[:, :, :, 1]

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
    # Remove the entire axis
    ax.axis('off')

    # plot the init frame
    p = [ax.plot([x[i, frame_idx]], [y[i, frame_idx]], c=colors[i], marker='o', markersize=6, linestyle='None', alpha=0.6) for i in range(n_agent)]

    fov = []
    for i in range(n_agent):
        pos = [x[i, frame_idx], y[i, frame_idx]]
        fov_x, fov_y = fov_xy(pos, yaw[i, frame_idx], fov_beta, fov_range)
        fov.append(ax.fill(np.concatenate(([pos[0]], fov_x, [pos[0]])), np.concatenate(([pos[1]], fov_y, [pos[1]])), color=colors[i], alpha=0.2))

    # Add the rectangles (boxes)
    boxes = []
    for i in range(n_agent):
        rect = plt.Rectangle((x[i][frame_idx] - bbox[0], y[i][frame_idx] - bbox[1]), 2*bbox[0], 2*bbox[1],
                             linewidth=1, edgecolor='blue', facecolor='none')
        ax.add_patch(rect)
        boxes.append(rect)

    if Trailing:
        trail = [Line2D([x[0][0]], [y[0][0]]) for i in range(n_agent)]
        for i in range(n_agent):
                trail[i], = ax.plot(x[i, :frame_idx+1], y[i, :frame_idx+1], c=colors[i], marker=None, linewidth=3, alpha=0.6)

    pred_idx = int(frame_idx//int(dt/Ts))
    if Estimation:
        # est_mean = [[Line2D([x[i][j]], [y[i][j]]) for j in range(n_agent-1)] for i in range(n_agent)]
        # for i in range(n_agent):
        #     for j in range(n_agent-1):
        #         est_mean[i][j], = ax.plot([est_x[i, j, pred_idx], est_y[i, j, pred_idx]], c=colors[i], marker=None, alpha=0.2)

        ellipses = []
        for i in range(n_agent):
            ellipses.append([])
            for j in range(n_agent-1):
                est_pos = estimate_mean[i, j, pred_idx, :2]  # [2, ]
                est_cov_arr = estimate_cov[i, j, pred_idx, :]  # [9, ]
                est_cov = np.array([est_cov_arr[0], est_cov_arr[1], est_cov_arr[3], est_cov_arr[4]]).reshape(2,2)
                ellipse_ij = plot_cov_ellipse_2d(est_pos, est_cov, edgecolor=colors[i], facecolor=colors[i], alpha=0.3, ax=ax)
                ellipses[i].append(ellipse_ij)

        p_near_ellipse = []
        for i in range(n_agent):
            p_near_ellipse.append([])
            for j in range(n_agent-1):
                try:
                    p_near_ellipse[i].append(ax.plot(p_near[i,j,pred_idx,0], p_near[i,j,pred_idx,1], c=colors[i], marker='x', markersize=5, linestyle='None', alpha=0.6))
                except:
                    print()

    if PredCurve:
        fov_preview_plot_step = list(range(0, pred_curve.shape[-2], preview_skip_step))+list([pred_curve.shape[-2]-1])
        iter_color = ['b', 'r', 'y', 'purple', 'pink', 'c', 'green', 'k']
        impc_iters = pred_curve.shape[-3]

        preds = []
        pred_fovs = []  # Store predicted FoV polygons
        for i in range(n_agent):
            iter_pred = []
            iter_pred_fov = []
            for impc_it in range(impc_iters):
                iter_pred_curve = pred_curve[:,:,impc_it,:,:]
                pred, = ax.plot(iter_pred_curve[i, pred_idx, :, 0], iter_pred_curve[i, pred_idx, :, 1], c=iter_color[impc_it], marker=None, alpha=0.8)
                iter_pred.append(pred)
                if enable_preview:
                    if impc_it == impc_iters-1 or impc_it == 0:
                        pred_fov = []
                        for k in fov_preview_plot_step:
                            pos = [iter_pred_curve[i, pred_idx, k, 0], iter_pred_curve[i, pred_idx, k, 1]]
                            fov_x, fov_y = fov_xy(pos, iter_pred_curve[i, pred_idx, k, 2], fov_beta, preview_fov_range)
                            pred_fov.append(ax.fill(np.concatenate(([pos[0]], fov_x, [pos[0]])),
                                                    np.concatenate(([pos[1]], fov_y, [pos[1]])),
                                                    color=iter_color[impc_it], alpha=0.1))  # Predicted FoV in red
                        iter_pred_fov.append(pred_fov)
            preds.append(iter_pred)
            pred_fovs.append(iter_pred_fov)


    if ShowGoals:
        for i in range(n_agent):
            goal = goals[i]
            circle = Circle((goal[0], goal[1]), goal_radius, alpha=0.2, linewidth=2)
            # Add the circle patch to the axes
            ax.add_patch(circle)

    plt.savefig(pre_save_name+"_frame"+str(frame_idx)+".pdf", bbox_inches='tight')
    # plt.show()


def animation2D_XYYaw(traj, estimate_mean, estimate_cov, p_near, dt, Ts, bbox, pred_curve=None, fov_beta=None, fov_range=None, Trailing=True, Estimation=True, PredCurve=True, colors="r", save_name= "./test.mp4"):
    # animation settings
    n_agent, total_frame, _ = traj.shape
    trailing_buf_size = 1000000
    dir_len = 0.1  # length for orientation tick

    # load traj
    states = traj  # [n, ts, dim]
    x = traj[:, :, 0]
    y = traj[:, :, 1]
    yaw = traj[:, :, 2]
    x_v = traj[:, :, 3]
    y_v = traj[:, :, 4]
    yaw_v = traj[:, :, 5]
    if Estimation:
        est_x = estimate_mean[:, :, :, 0]  # [n, n-1, ts]
        est_y = estimate_mean[:, :, :, 1]

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
    p = [ax.plot(x[0], y[0], c=colors[i], marker='o', markersize=6, linestyle='None', alpha=0.6) for i in range(n_agent)]

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
        fov.append(ax.fill(np.concatenate(([pos[0]], fov_x, [pos[0]])), np.concatenate(([pos[1]], fov_y, [pos[1]])), color=colors[i], alpha=0.2))

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
            trail[i], = ax.plot([x[i, 0]], [y[i, 0]], c=colors[i], marker=None, linewidth=3, alpha=0.4)

    if Estimation:
        est_mean = [[Line2D([x[i][j]], [y[i][j]]) for j in range(n_agent-1)] for i in range(n_agent)]
        for i in range(n_agent):
            for j in range(n_agent-1):
                est_mean[i][j], = ax.plot([est_x[i, j, 0], est_y[i, j, 0]], c=colors[i], marker=None, alpha=0.2)

        ellipses = []
        for i in range(n_agent):
            ellipses.append([])
            for j in range(n_agent-1):
                est_pos = estimate_mean[i, j, 0, :2]  # [2, ]
                est_cov_arr = estimate_cov[i, j, 0, :]  # [9, ]
                est_cov = np.array([est_cov_arr[0], est_cov_arr[1], est_cov_arr[3], est_cov_arr[4]]).reshape(2,2)
                ellipse_ij = plot_cov_ellipse_2d(est_pos, est_cov, edgecolor=colors[i], facecolor=colors[i], alpha=0.3, ax=ax)
                ellipses[i].append(ellipse_ij)

        p_near_ellipse = []
        for i in range(n_agent):
            p_near_ellipse.append([])
            for j in range(n_agent-1):
                p_near_ellipse[i].append(ax.plot(p_near[i,j,0,0], p_near[i,j,0,1], c=colors[i], marker='x', markersize=5, linestyle='None', alpha=0.6))

    if PredCurve:
        fov_preview_plot_step = list(range(0, pred_curve.shape[-2], preview_skip_step))+list([pred_curve.shape[-2]-1])
        iter_color = ['b', 'r', 'y', 'purple', 'pink', 'c', 'green', 'k']
        impc_iters = pred_curve.shape[-3]

        preds = []
        pred_fovs = []  # Store predicted FoV polygons
        for i in range(n_agent):
            iter_pred = []
            iter_pred_fov = []
            for impc_it in range(impc_iters):
                iter_pred_curve = pred_curve[:,:,impc_it,:,:]
                pred, = ax.plot([x[i, 0]], [y[i, 0]], c=iter_color[impc_it], marker=None, alpha=0.8)
                iter_pred.append(pred)
                if enable_preview:
                    if impc_it == impc_iters-1 or impc_it == 0:
                        pred_fov = []
                        for k in fov_preview_plot_step:
                            pos = [iter_pred_curve[i, 0, k, 0], iter_pred_curve[i, 0, k, 1]]
                            fov_x, fov_y = fov_xy(pos, iter_pred_curve[i, 0, k, 2], fov_beta, preview_fov_range)
                            pred_fov.append(ax.fill(np.concatenate(([pos[0]], fov_x, [pos[0]])),
                                                    np.concatenate(([pos[1]], fov_y, [pos[1]])),
                                                    color=iter_color[impc_it], alpha=0.1))  # Predicted FoV in red
                        iter_pred_fov.append(pred_fov)
            preds.append(iter_pred)
            pred_fovs.append(iter_pred_fov)


    # loop the frame
    def animate(ts):
        if ts % 100 == 0:
            print("Animating frame {0}/{1}".format(ts, total_frame))

        for i in range(n_agent):
            p_i = p[i][0]
            p_i.set_data([x[i, ts]], [y[i, ts]])
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
        if Estimation:
            for i in range(n_agent):
                for j in range(n_agent-1):
                    est_mean[i][j].set_data([est_x[i, j, :pred_index], est_y[i, j, :pred_index]])

            for i in range(n_agent):
                for j in range(n_agent-1):
                    est_pos = estimate_mean[i, j, pred_index, :2]  # [2, ]
                    est_cov_arr = estimate_cov[i, j, pred_index, :]  # [9, ]
                    est_cov = np.array([est_cov_arr[0], est_cov_arr[1], est_cov_arr[3], est_cov_arr[4]]).reshape(2,2)

                    _, ellipse_width, ellipse_height, ellipse_theta = compute_ellipse(est_pos, est_cov)
                    ellipses[i][j].center = (est_pos[0], est_pos[1])
                    ellipses[i][j].width = ellipse_width
                    ellipses[i][j].height = ellipse_height
                    ellipses[i][j].angle = ellipse_theta

            for i in range(n_agent):
                for j in range(n_agent-1):
                    p_near_ellipse[i][j][0].set_data([p_near[i,j,pred_index,0]], [p_near[i,j,pred_index,1]])


        if PredCurve:
            for i in range(n_agent):
                for impc_it in range(impc_iters):
                    iter_pred_curve = pred_curve[:,:,impc_it,:,:]
                    preds[i][impc_it].set_data(iter_pred_curve[i, pred_index, :, 0], iter_pred_curve[i, pred_index, :, 1])
                    if enable_preview:
                        for index_i, iter_index in enumerate([0]): # TODO fix this, for now is a temporary fix and the next line is the original
                        # for index_i, iter_index in enumerate([0, impc_iters-1]):
                            if impc_it == iter_index:
                                # Update predicted FoV
                                for index, k in enumerate(fov_preview_plot_step):
                                    pos = [iter_pred_curve[i, pred_index, k, 0], iter_pred_curve[i, pred_index, k, 1]]
                                    fov_x, fov_y = fov_xy(pos, iter_pred_curve[i, pred_index, k, 2], fov_beta, preview_fov_range)
                                    pred_fovs[i][index_i][index][0].set_xy(np.column_stack((np.concatenate(([pos[0]], fov_x, [pos[0]])),
                                                                                            np.concatenate(([pos[1]], fov_y, [pos[1]])))))

        return p[0][0],

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, frames=total_frame, interval=Ts*1e+3, blit=True)
    # anim.save(save_name, writer=animation.FFMpegWriter(fps=1/Ts))
    anim.save(save_name, writer=animation.PillowWriter(fps=1/Ts))
    return anim

def plot2D_XYYaw(traj, goals, goal_radius=1, obs_time=None, save_name="./test.jpg"):

    n_agent = traj.shape[0]
    assert n_agent==goals.shape[0]
    fig, ax = plt.subplots()
    # Set equal aspect ratio
    ax.set_aspect('equal')
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
        ax.plot(x, y, label='motion' + str(obs_index))

    for i in range(n_agent):
        goal = goals[i]
        circle = Circle((goal[0], goal[1]), goal_radius, alpha=0.2, linewidth=2)
        # Add the circle patch to the axes
        ax.add_patch(circle)

    plt.legend()
    plt.title('2D Trajectory')
    plt.savefig(save_name)
    # plt.show()

def derivatives_plot(traj, dt):
    n_agent = traj.shape[0]
    ts = np.arange(traj.shape[1])
    fig, ax = plt.subplots()
    # Set equal aspect ratio
    # ax.set_aspect('equal')
    ax.set_ylim([-6, 6])

    x_v = traj[:, :, 3]
    x_y = traj[:, :, 4]
    x_yaw = traj[:, :, 5]


    for robot_idx in range(n_agent):
        ax.plot(ts, x_v[robot_idx])
        pass
    plt.show()

if __name__ == "__main__":
    default_instance_type="circle"
    default_instance=default_instance_type+"4"
    default_fov = 120
    exp_idx=0
    parser = argparse.ArgumentParser(
        description="argparse to read the config, states and output filenames"
    )
    # TODO: properly construct the inputs
    parser.add_argument("-c", "--config_filename", type=str, default="/usr/src/mpc-cbf/workspace/experiments/config/circle/circle2_config.json", help="path to config json file")
    parser.add_argument("-s", "--states_filename", type=str, default="/usr/src/mpc-cbf/workspace/experiments/results/states.json", help="path to simulation state json file")
    # Need to install ffmpeg to save the video as mp4 (use pillow to save as gif)
    parser.add_argument("-ov", "--output_video", type=str, default="../../results"+default_instance+".gif", help="path to simulation animation file")
    parser.add_argument("-of", "--output_figure", type=str, default="../../results"+default_instance, help="path to simulation figure file")
    parser.add_argument("-f", "--fov", type=int, default=default_fov, help="fov of simulation")
    args = parser.parse_args()

    config_json = args.config_filename
    states_json = load_states(args.states_filename)
    output_video = args.output_video
    output_figure = args.output_figure
    
    # make directory if not exist
    if not os.path.exists(os.path.dirname(output_video)):
        os.makedirs(os.path.dirname(output_video))
    if not os.path.exists(os.path.dirname(output_figure)):
        os.makedirs(os.path.dirname(output_figure))

    num_robots = len(states_json["robots"])
    traj = np.array([states_json["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
    # if "baseline" not in args.states_filename:
    #     traj = traj[:,::10,:]
    neighbor_ids = []
    for i in range(num_robots):
        neighbor_ids.append([])
        for j in range(num_robots):
            if i!=j:
                neighbor_ids[i].append(j)

    estimate_mean = []
    estimate_cov = []
    p_near = []
    for i in range(num_robots):
        i_neighbor_ids = neighbor_ids[i]
        estimate_mean.append([])
        estimate_cov.append([])
        p_near.append([])
        for j in i_neighbor_ids:
            estimate_mean[i].append(states_json["robots"][str(i)]["estimates_mean"][str(j)])
            estimate_cov[i].append(states_json["robots"][str(i)]["estimates_cov"][str(j)])
            p_near[i].append(states_json["robots"][str(i)]["p_near_ellipse"][str(j)])

    estimate_mean = np.array(estimate_mean)  # [n_robot, n_robot-1, ts, dim]
    estimate_cov = np.array(estimate_cov)  # [n_robot, n_robot-1, ts, dimxdim]
    p_near = np.array(p_near)  # [n_robot, n_robot-1, ts, 2]

    # estimate_mean = None
    # estimate_cov = None
    # p_near = None
    # pred_curve = None

    # estimate_mean = [[states_json["robots"][str(_)]["estimates_mean"][neighbor_id] for neighbor_id in neighbor_ids[_]] for _ in range(num_robots)]]  # [n_robot, n_robot-1, ts, dim]
    dt = states_json["dt"]
    Ts = states_json["Ts"]
    # if "baseline" not in args.states_filename:
    #     Ts = Ts*10
    bbox = load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2]
    goals = np.array(load_states(config_json)["tasks"]["sf"])
    goal_radius = 1.6

    # pred_curve = [states_json["robots"][str(_)]["pred_curve"] for _ in range(num_robots)]  # [n, h_samples, impc_iter(dynamic), horizon, dim]
    # impc_iter_max = max(
    #     max(len(seq) for seq in sublist if seq is not None)  # Ignore None values
    #     for sublist in pred_curve
    # )
    # d1 = len(pred_curve)
    # d2 = len(pred_curve[0])
    # d3 = impc_iter_max
    # for i in range(len(pred_curve[0])):
    #     if len(pred_curve[0][i]) != 0:
    #         d4 = len(pred_curve[0][i][0])
    #         d5 = len(pred_curve[0][i][0][0])
    #         break
    # # Create a numpy array filled with zeros (or another padding value)
    # padded_array = np.zeros((d1, d2, d3, d4, d5)) * np.nan
    # # Fill in the original data
    # for i in range(d1):
    #     for j in range(d2):
    #         seq = pred_curve[i][j]  # Access the varying `x` dimension
    #         if seq is not None and len(seq) > 0:  # Ensure it's not None
    #             seq_array = np.asarray(seq)  # Convert to NumPy array
    #             padded_array[i, j, :len(seq), :, :] = seq  # Copy actual data
    # pred_curve = padded_array
    # # # Fill the array with data
    # # for i, seq in enumerate(data_list):
    # # padded_array[i, :len(seq), :, :] = seq  # Copy data
    # # pred_curve = [states_json["robots"][str(_)]["pred_curve"] for _ in range(num_robots)]  # [n, h_samples, impc_iter(dynamic), horizon, dim]
    # pred_curve = np.array([states_json["robots"][str(_)]["pred_curve"] for _ in range(num_robots)])  # [n, h_samples, impc_iter, horizon, dim]
    pred_curve = None
    FoV_beta = args.fov * np.pi/180
    FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]

    colors = generate_rgb_colors(num_robots)

    # goal_radius = 1.6
    # success = instance_success(traj, goals, goal_radius, collision_shape)
    # print("success: ", success)
    derivatives_plot(traj, Ts)
    plot2D_XYYaw(traj, goals, save_name=output_figure)
    _, total_frame, _ = traj.shape
    frame_gap = 20
    render_frame_indices = range(0, total_frame, frame_gap)
    # for frame_idx in render_frame_indices:
    #     snapshots2D_XYYaw(traj, goals, goal_radius, estimate_mean, estimate_cov, p_near, dt, Ts, bbox, frame_idx, pred_curve, FoV_beta, FoV_range, PredCurve=True, colors=colors, pre_save_name=output_figure)
    animation2D_XYYaw(traj, estimate_mean, estimate_cov, p_near, dt, Ts, bbox, pred_curve, FoV_beta, FoV_range, Estimation=True, PredCurve=False, colors=colors, save_name=output_video)