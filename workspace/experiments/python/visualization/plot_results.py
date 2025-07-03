import argparse
import json
import os

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def generate_colors(n):
    hsv = np.linspace(0, 1, n + 1)[:-1]
    return [plt.cm.hsv(h) for h in hsv]

def plot_connectivity(ax, positions, max_dist, colors, label_prefix=''):
    """
    在 ax 坐标轴上绘制一次性的连通性：在 positions 中分别画出机器人位置，
    并对两两之间距离小于等于 max_dist 的机器人连线。
    positions: N×3 的数组 (x, y, θ)，这里只用前两列做平面坐标。
    colors: 长度为 N 的颜色列表，每个机器人一种颜色。
    """
    n = len(positions)
    for i in range(n):
        x_i, y_i = positions[i][:2]
        ax.plot(x_i, y_i, 'o', color=colors[i], label=f"{label_prefix}{i}")
        for j in range(i + 1, n):
            x_j, y_j = positions[j][:2]
            if np.hypot(x_j - x_i, y_j - y_i) <= max_dist:
                ax.plot([x_i, x_j], [y_i, y_j], '-', color='gray', linewidth=1)

def plot_trajectory(ax, traj, colors):
    """
    在 ax 上画出所有机器人的全程轨迹（静态），不会用于动画更新时调用。
    traj: 形状 (N, T, 3) 的数组，每辆机器人 T 帧的 (x, y, θ)。
    """
    n = traj.shape[0]
    for i in range(n):
        ax.plot(traj[i, :, 0], traj[i, :, 1], '-', color=colors[i], alpha=0.7)

def main():
    # === 1. 设定文件路径 ===
    parser = argparse.ArgumentParser(description="Plot formation connectivity and trajectories")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--states", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument(
        "--create_anim",
        action="store_true",
        help="Show animation in a window (requires GUI backend)",
    )

    args = parser.parse_args()
    config_file = args.config
    states_file = args.states

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    output_static = os.path.join(output_dir, os.path.basename(config_file).replace('.json', '.png'))
    output_anim = os.path.join(output_dir, os.path.basename(config_file).replace('.json', '.mp4'))

    # config_file = "/usr/src/mpc-cbf/workspace/experiments/config/circle/circle3_config.json"
    # states_file = "/usr/src/mpc-cbf/workspace/experiments/results/formation/states.json"
    # output_static = "/usr/src/mpc-cbf/workspace/experiments/results/formation/viz/plot_circle3.png"
    # output_anim = "/usr/src/mpc-cbf/workspace/experiments/results/formation/viz/anim_circle3.mp4"

    # === 2. 读取配置和状态数据 ===
    cfg = load_json(config_file)
    st = load_json(states_file)

    # 初始和最终的机器人数组
    sf = np.array(cfg['tasks']['sf'])
    so = np.array(cfg['tasks']['so'])
    max_dist = cfg['cbf_params']['d_max']

    # robots 数据结构：假设 st['robots'][key]['states'] 是一个帧序列
    robots = st['robots']
    keys = sorted(robots.keys(), key=int)  # 假设机器人 ID 是字符串数字，用 int 排序
    # 将每个机器人各帧状态组织成 (N, T, 3) 的数组
    traj = np.array([robots[k]['states'] for k in keys], dtype=float)  # 例如形状 (N, T, 3)
    N, T, _ = traj.shape

    colors = generate_colors(N)

    # === 3. 先画静态子图 ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # --- 3.1 初始连通性（Static） ---
    axes[0].set_title('Initial Connectivity')
    plot_connectivity(axes[0], so, max_dist, colors)
    axes[0].axis('equal')
    axes[0].grid(True)

    # --- 3.2 最终连通性（Static） ---
    axes[1].set_title('Final Connectivity')
    plot_connectivity(axes[1], sf, max_dist, colors)
    axes[1].axis('equal')
    axes[1].grid(True)

    # --- 3.3 轨迹 & （动画时会在这里叠加动态连通性） ---
    axes[2].set_title('Trajectories & Connectivity Over Time')
    plot_trajectory(axes[2], traj, colors)
    axes[2].axis('equal')
    axes[2].grid(True)

    # 保存静态图
    os.makedirs(os.path.dirname(output_static) or '.', exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_static)
    print(f"Static plot saved to {output_static}")

    # === 4. 准备动画：在第三个子图 axes[2] 上，动态绘制每帧的连通性 ===
    # 为了在每帧更新时能删除前一帧绘制的连通性线、散点等，需要先保存第三个子图 axes[2] 的背景
    ax_anim = axes[2]

    # 我们会在动画里，每帧清除上一次绘制的 “连通性标记”（scatter + lines），
    # 但保留轨迹线（因为轨迹用 plot 绘制，多条线不必反复重画）。
    # 最简单的做法是为每帧重新在 axes[2] 上画机器人当前位置和连通性。

    # 定义一个存放每一帧绘图对象的列表，用以 update 时清理
    artists = []

    def init_frame():
        """
        初始化函数：在动画开头就画出全程轨迹（已经在静态图绘制过），
        这里只负责“清空”第三个子图上可动的艺术家对象。
        """
        # 把前一次留下的所有 artists 都 remove 掉
        nonlocal artists
        for art in artists:
            try:
                art.remove()
            except Exception:
                pass
        artists = []
        return []

    def update_frame(frame_idx):
        """
        更新函数：每次动画推进到第 frame_idx 帧时调用
        frame_idx 范围 [0, T-1]：
        - 在 (x, y) 处画 scatter 点（可以用小圆圈表示当前帧机器人位置）
        - 检查两两距离，画连通性连线（灰色线）
        """
        nonlocal artists
        # 先把上一次画的动态对象都删掉
        for art in artists:
            try:
                art.remove()
            except Exception:
                pass
        artists = []

        # 取出当前帧每辆车的 (x, y)
        curr_positions = traj[:, frame_idx, :2]  # 形状 (N, 2)
        # 先用 scatter 画机器人当前位置
        scat = ax_anim.scatter(
            curr_positions[:, 0],
            curr_positions[:, 1],
            c=colors,
            s=50,
            edgecolors='k',
            zorder=3
        )
        artists.append(scat)

        # 再两两判断距离，画连通性线段
        for i in range(N):
            x_i, y_i = curr_positions[i]
            for j in range(i + 1, N):
                x_j, y_j = curr_positions[j]
                if np.hypot(x_j - x_i, y_j - y_i) <= max_dist:
                    ln, = ax_anim.plot([x_i, x_j], [y_i, y_j], '-', color='gray', lw=1, zorder=2)
                    artists.append(ln)

        return artists

    # 构造动画对象：共 T 帧，interval=200ms

    if args.create_anim:
        anim = animation.FuncAnimation(
            fig,
            update_frame,
            frames=T,
            init_func=init_frame,
            blit=False,  # blit=True 在某些环境下需额外处理 background；如果卡，可改为 False
            interval=200,  # 每帧间隔毫秒数，可根据需要调整
            repeat=False,
        )

        # === 5. 保存动画 ===
        # 如果要保存为 mp4，需要系统上安装了 ffmpeg 或者 avconv
        os.makedirs(os.path.dirname(output_anim) or ".", exist_ok=True)
        try:
            # 保存为 MP4，帧率为 5 帧/秒 (fps=1000/interval ≈ 5)
            anim.save(output_anim, writer="ffmpeg", fps=5)
            print(f"Animation saved to {output_anim}")
        except Exception as e:
            print("Failed to save as MP4: ", e)
            # 如果想保存为 GIF，可以尝试下面的方式（需 imagemagick 支持）：
            gif_path = output_anim.replace(".mp4", ".gif")
            try:
                anim.save(gif_path, writer="imagemagick", fps=5)
                print(f"Animation also saved to {gif_path}")
            except Exception as e2:
                print("Also failed to save as GIF:", e2)

if __name__ == '__main__':
    main()
