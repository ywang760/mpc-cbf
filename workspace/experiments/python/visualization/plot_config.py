# Given a configuration file, extract that starting and goal positions

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot configuration from JSON file")
    parser.add_argument("config_file", type=str, help="Path to the configuration JSON file")
    args = parser.parse_args()

    # Load the configuration file
    with open(args.config_file, 'r') as f:
        config = json.load(f)

    # Extract start and goal positions
    start_positions = np.array(config['tasks']['so'])
    goal_positions = np.array(config['tasks']['sf'])

    # Create two plots side by side
    num_robots = len(start_positions)

    # Generate unique colors for each robot
    colors = plt.cm.viridis(np.linspace(0, 1, num_robots))

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # Plot start positions with unique colors per robot
    for i in range(num_robots):
        axs[0].scatter(start_positions[i, 0], start_positions[i, 1], 
                      color=colors[i], s=60, label=f'Robot {i+1}')
    axs[0].set_title('Start Positions')
    axs[0].set_xlabel('X Coordinate')
    axs[0].set_ylabel('Y Coordinate')
    axs[0].grid(True)
    axs[0].legend()

    # Plot goal positions with the same colors per robot
    for i in range(num_robots):
        axs[1].scatter(goal_positions[i, 0], goal_positions[i, 1], 
                      color=colors[i], s=60, label=f'Robot {i+1}')
    axs[1].set_title('Goal Positions')
    axs[1].set_xlabel('X Coordinate')
    axs[1].set_ylabel('Y Coordinate')
    axs[1].grid(True)
    axs[1].legend()

    # Set the same axis limits for both plots to enable direct comparison
    all_positions = np.vstack([start_positions, goal_positions])

    x_min, y_min = config["physical_limits"]["p_min"]
    x_max, y_max = config["physical_limits"]["p_max"]
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.1

    for ax in axs:
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
        ax.set_aspect('equal', adjustable='box')

    plt.suptitle('Start and Goal Positions from Configuration File')

    # /usr/src/mpc-cbf/workspace/experiments/c^Cfig/formation/robots2_circle.
    save_dir = os.path.dirname(args.config_file)
    save_dir = os.path.join(save_dir, 'viz')
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, os.path.basename(args.config_file).replace('.json', '.png'))
    print(f"Saving plot to {save_path}")
    plt.savefig(save_path)
