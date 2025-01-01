import numpy as np
import matplotlib.pyplot as plt
import json

def generate_points_on_circle(num_points, radius, angle_bias):
    """
    Generates (num_points) points uniformly at random on the circumference
    of a circle with a given (radius).

    Returns:
        x_coords (ndarray): Array of x-coordinates.
        y_coords (ndarray): Array of y-coordinates.
    """
    # Incremental angles from [angle_bias, angle_bias+2pi]
    thetas = np.linspace(angle_bias, angle_bias+2*np.pi, num=num_points, endpoint=False)

    # Convert polar to Cartesian
    x_coords = radius * np.cos(thetas)
    y_coords = radius * np.sin(thetas)

    return x_coords, y_coords

def compute_yaw(x, y, circle_center):
    yaws = np.arctan2(-y - circle_center[1], -x - circle_center[0])
    return yaws

def fov_xy(x, y, radian, fov_beta, fov_range, num_points=100):
    pos_x = x
    pos_y = y
    angles = np.linspace(-fov_beta/2 + radian, fov_beta/2 + radian, num_points)
    # Calculate the points on the edges of the FOV
    x = fov_range * np.cos(angles) + pos_x
    y = fov_range * np.sin(angles) + pos_y
    return x,y

def plot_fov(x, y, fov_beta, fov_range, yaw, c="r", ax=None):
    fov_x, fov_y = fov_xy(x, y, yaw, fov_beta, fov_range)
    ax.fill(np.concatenate(([x], fov_x, [x])),
            np.concatenate(([y], fov_y, [y])),
            color=c, alpha=0.3)

def plot_position(x, y, idx, c="r", ax=None):
    ax.scatter(x, y, c=c)
    ax.text(x, y, str(idx), fontsize=10)

if __name__ == '__main__':
    # mpc params
    mpc_params = {
        "h": 0.1,
        "Ts": 0.01,
        "k_hor": 16,
        "mpc_tuning": {
            "w_pos_err": 10,
            "w_u_eff": 1,
            "spd_f": 3
        },
        "physical_limits": {
            "p_min": [-10, -10],
            "p_max": [10, 10],
            "v_min": [-0.5, -0.5, -0.5],
            "v_max": [0.5, 0.5, 0.5],
            "a_min": [-10.0, -10.0, -10.0],
            "a_max": [10.0, 10.0, 10.0],
            "pos_std": 0.001,
            "vel_std": 0.01
        }
    }

    # bezier_params
    bezier_params = {
        "num_pieces": 3,
        "num_control_points": 4,
        "piece_max_parameter": 0.5
    }

    # fov cbf params
    fov_cbf_params = {
        "beta": 120,
        "Rs": 10
    }

    # robot params
    robot_params = {
        "collision_shape": {
            "aligned_box": [0.3, 0.3, 0]
        }
    }

    # position config
    circle_radius = 4
    circle_center = np.array([0, 0])
    num_robots = 8
    angle_bias = 0.1

    start_xs, start_ys = generate_points_on_circle(num_robots, circle_radius, angle_bias)  # [num_robots,], [num_robots,]
    goal_xs, goal_ys = generate_points_on_circle(num_robots, circle_radius, np.pi-angle_bias)  # [num_robots,], [num_robots,]

    start_yaws = compute_yaw(start_xs, start_ys, circle_center)
    goal_yaws = compute_yaw(goal_xs, goal_ys, circle_center)
    # print("start_yaws: ", start_yaws)
    # print("goal_yaws: ", goal_yaws)

    so = []
    sf = []
    for i in range(num_robots):
        so.append([start_xs[i], start_ys[i], start_yaws[i]])
        sf.append([goal_xs[i], goal_ys[i], goal_yaws[i]])

    # visualize the position and fov
    fov_beta = 120 * np.pi / 180
    fov_range = 1
    fig, ax = plt.subplots()
    # Set equal aspect ratio
    ax.set_aspect('equal')
    for i in range(num_robots):
        plot_position(start_xs[i], start_ys[i], i, "r", ax)
        plot_fov(start_xs[i], start_ys[i], fov_beta, fov_range, start_yaws[i], "r", ax)
        plot_position(goal_xs[i], goal_ys[i], i, "b", ax)
        plot_fov(goal_xs[i], goal_ys[i], fov_beta, fov_range, goal_yaws[i], "b", ax)

    # plt.plot()
    plt.show()

    # generate the json config file
    data = {
        "mpc_params": mpc_params,
        "bezier_params": bezier_params,
        "fov_cbf_params": fov_cbf_params,
        "tasks": {
            "so": so,
            "sf": sf
        },
        "robot_params": robot_params
    }
    # save to json
    with open("circle%d_config.json"%(num_robots), "w") as file:
        json.dump(data, file, indent=4)