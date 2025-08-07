import argparse
import json
import os
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

# Suppress matplotlib animation warnings
warnings.filterwarnings('ignore', message='MovieWriter.*unavailable')


def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def generate_colors(n):
    hsv = np.linspace(0, 1, n + 1)[:-1]
    return [plt.cm.hsv(h) for h in hsv]

def plot_connectivity(ax, positions, max_dist, colors, radius, label_prefix=''):
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
        circle = plt.Circle((x_i, y_i), radius, color=colors[i], alpha=0.2)
        ax.add_patch(circle)
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
    parser = argparse.ArgumentParser(description="Plot formation connectivity and trajectories")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--states", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument(
        "--create_anim",
        action="store_true",
    )
    parser.add_argument(
        "--anim_format",
        type=str,
        choices=["mp4", "gif"],
        default="gif",
        help="Format for the animation output (default: gif)",
    )

    args = parser.parse_args()
    config_file = args.config
    states_file = args.states

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    output_dir_expanded = os.path.join(
        output_dir, os.path.dirname(config_file).split("/")[-1]
    )
    output_static = os.path.join(
        output_dir_expanded, os.path.basename(config_file).replace(".json", ".png")
    )
    output_anim = os.path.join(
        output_dir_expanded, os.path.basename(config_file).replace(".json", ".mp4")
    )

    # === 2. 读取配置和状态数据 ===
    cfg = load_json(config_file)
    robot_radius = cfg["robot_params"]["collision_shape"]["radius"]
    # try to load states.json, if it doesn't exist, it will raise an error
    HAS_RESULT = False
    if os.path.exists(states_file):
        st = load_json(states_file)
        HAS_RESULT = True

    # 初始和最终的机器人数组
    sf = np.array(cfg['tasks']['sf'])
    so = np.array(cfg['tasks']['so'])
    max_dist = cfg['cbf_params']['d_max']
    ts = cfg["mpc_params"]["h"]
    N = len(sf)
    

    # robots 数据结构：假设 st['robots'][key]['states'] 是一个帧序列
    T = 0  # Initialize T to avoid unbound variable issues
    if HAS_RESULT:
        robots = st["robots"]
        keys = sorted(robots.keys(), key=int)  # 假设机器人 ID 是字符串数字，用 int 排序
        # 将每个机器人各帧状态组织成 (N, T, 3) 的数组
        traj = np.array(
            [robots[k]["states"] for k in keys], dtype=float
        )  # 例如形状 (N, T, 3)
        _, T, _ = traj.shape

    colors = generate_colors(N)

    # === 3. 创建子图 ===
    if HAS_RESULT:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # --- 3.1 初始连通性（Static） ---
    axes[0].set_title('Initial Connectivity')
    plot_connectivity(axes[0], so, max_dist, colors, robot_radius)
    axes[0].grid(True)

    # --- 3.2 最终连通性（Static） ---
    axes[1].set_title('Final Connectivity')
    plot_connectivity(axes[1], sf, max_dist, colors, robot_radius)
    axes[1].grid(True)

    # --- 3.3 轨迹（仅当有结果时） ---
    if HAS_RESULT:
        axes[2].set_title("Trajectories & Connectivity Over Time")
        plot_trajectory(axes[2], traj, colors)
        axes[2].grid(True)

    # 保存静态图
    # Set x and y limits based on physical limits
    x_min, y_min = cfg["physical_limits"]["p_min"]
    x_max, y_max = cfg["physical_limits"]["p_max"]
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.1
    for ax in axes:
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
    os.makedirs(os.path.dirname(output_static) or '.', exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_static)
    print(f"Static plot saved to {output_static}")

    if args.create_anim and HAS_RESULT:
        
        ax_anim = axes[2]
        artists = []
        
        def init_frame():
            nonlocal artists
            for art in artists:
                try:
                    art.remove()
                except Exception:
                    pass
            artists = []
            return []

        def update_frame(frame_idx):
            nonlocal artists
            for art in artists:
                try:
                    art.remove()
                except Exception:
                    pass
            artists = []

            # 取出当前帧每辆车的 (x, y)
            curr_positions = traj[:, frame_idx, :2]  # 形状 (N, 2)

            # 1. scatter画出当前位置
            scat = ax_anim.scatter(
                curr_positions[:, 0],
                curr_positions[:, 1],
                c=colors,
                s=50,
                edgecolors='k',
                zorder=3
            )
            artists.append(scat)

            # 2. 每个机器人画圆
            for i in range(N):
                x_i, y_i = curr_positions[i]
                circle = plt.Circle((x_i, y_i), robot_radius, color=colors[i], alpha=0.2, zorder=1)
                ax_anim.add_patch(circle)
                artists.append(circle)

            # 3. 两两之间画连线
            for i in range(N):
                x_i, y_i = curr_positions[i]
                for j in range(i + 1, N):
                    x_j, y_j = curr_positions[j]
                    if np.hypot(x_j - x_i, y_j - y_i) <= max_dist:
                        ln, = ax_anim.plot([x_i, x_j], [y_i, y_j], '-', color='gray', lw=1, zorder=2)
                        artists.append(ln)

            return artists


        # Downsample frames if T is too large for smooth animation
        MAX_FRAMES = 200  # Maximum frames for reasonable animation
        if T > MAX_FRAMES:
            # Calculate downsampling factor
            downsample_factor = max(1, T // MAX_FRAMES)
            frame_indices = np.arange(0, T, downsample_factor)
            actual_frames = len(frame_indices)
            actual_interval = 1000 * ts * downsample_factor  # Adjust interval for downsampling
            print(f"Downsampling animation: {T} frames -> {actual_frames} frames (factor: {downsample_factor})")
            
            def update_frame_downsampled(frame_idx):
                # Use the downsampled frame index
                actual_frame_idx = frame_indices[frame_idx]
                return update_frame(actual_frame_idx)
            
            anim = animation.FuncAnimation(
                fig,
                update_frame_downsampled,
                frames=actual_frames,
                init_func=init_frame,
                blit=False,
                interval=actual_interval,
                repeat=False,
            )
        else:
            anim = animation.FuncAnimation(
                fig,
                update_frame,
                frames=T,
                init_func=init_frame,
                blit=False,
                interval=1000 * ts,
                repeat=False,
            )

        # === 5. 保存动画 ===
        os.makedirs(os.path.dirname(output_anim) or ".", exist_ok=True)

        # Set output path and writer based on format
        if args.anim_format == "mp4":
            output_path = output_anim
            writer = "ffmpeg"
        else:  # gif
            output_path = output_anim.replace(".mp4", ".gif")
            writer = "Pillow"

        try:
            anim.save(output_path, writer=writer, fps=1 / ts)
            print(f"Animation saved to {output_path}")
        except Exception as e:
            print(f"Failed to save as {args.anim_format.upper()}: {e}")
    else:
        print("Animation creation skipped or no results available.")


if __name__ == '__main__':
    main()
