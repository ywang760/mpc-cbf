import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os

def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)


def generate_colors(n):
    hsv = np.linspace(0, 1, n + 1)[:-1]
    return [plt.cm.hsv(h) for h in hsv]


def compute_ellipse_params(cov2d, nstd=2):
    # eigenvalues and eigenvectors
    vals, vecs = np.linalg.eigh(cov2d)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2 * nstd * np.sqrt(vals)
    return width, height, theta


def plot_connectivity(ax, positions, max_dist, colors, label_prefix=''):
    n = len(positions)
    for i in range(n):
        x_i, y_i = positions[i][:2]
        ax.plot(x_i, y_i, 'o', color=colors[i], label=f"{label_prefix}{i}")
        for j in range(i + 1, n):
            x_j, y_j = positions[j][:2]
            if np.hypot(x_j - x_i, y_j - y_i) <= max_dist:
                ax.plot([x_i, x_j], [y_i, y_j], '-', color='gray', linewidth=1)


def plot_trajectory(ax, traj, colors):
    n = traj.shape[0]
    for i in range(n):
        ax.plot(traj[i, :, 0], traj[i, :, 1], '-', color=colors[i], alpha=0.7)


def plot_cov_ellipses(ax, states, covs, colors, step=10):
    n = len(states)
    for i in range(n):
        for j in range(len(states[i])):
            if j % step != 0:
                continue
            pos = states[i][j][:2]
            cov9 = covs[i][j]
            cov2 = np.array([[cov9[0], cov9[1]], [cov9[3], cov9[4]]])
            w, h, ang = compute_ellipse_params(cov2)
            ellipse = Ellipse(xy=pos, width=w, height=h, angle=ang,
                              edgecolor=colors[i], facecolor='none', alpha=0.3)
            ax.add_patch(ellipse)


def main():
    config_file = "/usr/src/mpc-cbf/workspace/experiments/config/circle/circle3_config.json"
    states_file = "/usr/src/mpc-cbf/workspace/experiments/results/states.json"
    output_file = "/usr/src/mpc-cbf/workspace/experiments/results"

    cfg = load_json(config_file)
    st = load_json(states_file)
    sf = np.array(cfg['tasks']['sf'])

    so = np.array(cfg['tasks']['so'])
    sf = np.array(cfg['tasks']['sf'])
    max_dist = cfg['connectivity_params']['max_distance']

    robots = st['robots']
    keys = sorted(robots.keys(), key=int)
    traj = np.array([robots[k]['states'] for k in keys], dtype=float)
    # covs = [robots[k]['estimates_cov'][keys[j]] for j,k in enumerate(keys)]
    # covs = [[np.array(c) for c in robots[k]['estimates_cov'][k2]]
    #         for k in keys for k2 in [keys[(int(k)+1)%len(keys)]]]  # simplified

    colors = generate_colors(len(keys))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # initial connectivity
    axes[0].set_title('Initial Connectivity')
    plot_connectivity(axes[0], so, max_dist, colors)
    axes[0].axis('equal')

    # final connectivity
    axes[1].set_title('Final Connectivity')
    plot_connectivity(axes[1], sf, max_dist, colors)
    axes[1].axis('equal')

    # trajectories and covariances
    axes[2].set_title('Trajectories & Covariance')
    plot_trajectory(axes[2], traj, colors)
    # you can add ellipses with plot_cov_ellipses if desired
    axes[2].axis('equal')

    for ax in axes:
        ax.grid(True)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    plt.savefig(output_file)
    print(f"Saved plot to {output_file}")

if __name__ == '__main__':
    main()
