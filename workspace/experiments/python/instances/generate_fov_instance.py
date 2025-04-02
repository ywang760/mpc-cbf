import numpy as np
import matplotlib.pyplot as plt
import json
import argparse
import os
import sys
from utils import *


def fov_xy(x, y, radian, fov_beta, fov_range, num_points=100):
    pos_x = x
    pos_y = y
    angles = np.linspace(-fov_beta / 2 + radian, fov_beta / 2 + radian, num_points)
    # Calculate the points on the edges of the FOV
    x = fov_range * np.cos(angles) + pos_x
    y = fov_range * np.sin(angles) + pos_y
    return x, y

def plot_fov(x, y, fov_beta, fov_range, yaw, c="r", ax=None):
    fov_x, fov_y = fov_xy(x, y, yaw, fov_beta, fov_range)
    ax.fill(
        np.concatenate(([x], fov_x, [x])),
        np.concatenate(([y], fov_y, [y])),
        color=c,
        alpha=0.3,
    )

if __name__ == "__main__":
    # take argument in
    parser = argparse.ArgumentParser(
        description="argparse to read the number of robots"
    )
    parser.add_argument(
        "-n",
        "--num_robots",
        type=int,
        required=True,
        help="number of robots in simulation",
    )
    parser.add_argument(
        "-ins",
        "--instance_type",
        type=str,
        default="circle",
        choices=["circle", "formation"],
        help="instance type (circle or formation)",
    )
    parser.add_argument("-f", "--fov", type=int, default=120, help="degree of fov")
    args = parser.parse_args()

    if args.num_robots <= 0:
        parser.error("Number of robots must be positive")
    if not 0 <= args.fov <= 360:
        parser.error("FOV must be between 0 and 360 degrees")
        
    # Load default parameters
    default_config_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "../../config/params/default_params.json"
    )
    with open(default_config_path, "r") as file:
        default_config = json.load(file)    
        
    # fov config
    instance_type = args.instance_type
    if instance_type == "circle":
        w_pos_err = 10
    elif instance_type == "formation":
        w_pos_err = 300

    # mpc params
    mpc_params = default_config["mpc_params"]
    mpc_params["mpc_tuning"]["w_pos_err"] = w_pos_err

    # fov cbf params
    fov_cbf_params = {"beta": args.fov, "Rs": 20}

    so = []
    sf = []
    num_robots = args.num_robots
    if instance_type == "circle":
        # position config
        circle_radius = 4
        circle_center = np.array([0, 0])
        angle_bias = 0.1

        start_xs, start_ys = generate_points_on_circle(
            num_robots, circle_radius, angle_bias
        )  # [num_robots,], [num_robots,]
        goal_xs, goal_ys = generate_points_on_circle(
            num_robots, circle_radius, np.pi - angle_bias
        )  # [num_robots,], [num_robots,]

        start_yaws = compute_yaw(start_xs, start_ys, circle_center)
        goal_yaws = compute_yaw(goal_xs, goal_ys, circle_center)

        for i in range(num_robots):
            so.append([start_xs[i], start_ys[i], start_yaws[i]])
            sf.append([goal_xs[i], goal_ys[i], goal_yaws[i]])

    elif instance_type == "formation":
        center_offset = 6
        formation_distance_x = 1
        formation_distance_y = 1
        n_row = 2
        start_center = np.array([-center_offset, 0])
        goal_center = np.array([center_offset, 0])

        start_xs, start_ys = generate_points_on_formation(
            num_robots, n_row, formation_distance_x, formation_distance_y, start_center
        )
        goal_xs, goal_ys = generate_points_on_formation(
            num_robots, n_row, formation_distance_x, formation_distance_y, goal_center
        )

        start_yaws = np.repeat(0.0, num_robots)
        goal_yaws = np.repeat(0.0, num_robots)

        for i in range(num_robots):
            so.append([start_xs[i], start_ys[i], start_yaws[i]])
            sf.append([goal_xs[i], goal_ys[i], goal_yaws[i]])

    # visualize the position and fov
    fov_beta = 120 * np.pi / 180
    fov_range = 1
    fig, ax = plt.subplots()
    # Set equal aspect ratio
    ax.set_aspect("equal")
    for i in range(num_robots):
        plot_position(start_xs[i], start_ys[i], i, "r", ax)
        plot_fov(start_xs[i], start_ys[i], fov_beta, fov_range, start_yaws[i], "r", ax)
        plot_position(goal_xs[i], goal_ys[i], i, "b", ax)
        plot_fov(goal_xs[i], goal_ys[i], fov_beta, fov_range, goal_yaws[i], "b", ax)
    plt.show()

    # generate the json config file
    data = {
        "mpc_params": mpc_params,
        "bezier_params": default_config["bezier_params"],
        "fov_cbf_params": fov_cbf_params,
        "robot_params": default_config["robot_params"],
        "tasks": {
            "so": so,
            "sf": sf,
        }
    }
    
    # save to json
    output_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), f"../../config/{instance_type}_instances"
    )
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(os.path.join(output_dir, f"{instance_type}{num_robots}_config.json"), "w") as file:
        json.dump(data, file, indent=4)
