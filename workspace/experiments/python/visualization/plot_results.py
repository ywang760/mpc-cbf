import argparse
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 触发 3D 投影支持


def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def generate_colors(n):
    hsv = np.linspace(0, 1, n + 1)[:-1]
    return [plt.cm.hsv(h) for h in hsv]

def plot_connectivity(ax, positions, max_dist, colors, radius, label_prefix='', is3d=False):
    """
    绘制连通性：
    - is3d=False: 用 (x,y)
    - is3d=True : 用 (x,y,z) 且距离用 3D 范数
    """
    n = len(positions)
    for i in range(n):
        if is3d:
            x_i, y_i, z_i = positions[i][:3]
            # 点
            ax.scatter([x_i], [y_i], [z_i], color=colors[i], label=f"{label_prefix}{i}")
            # z=z_i 平面的“半径环”（3D 无法直接加 2D Circle patch）
            t = np.linspace(0, 2*np.pi, 80)
            ring_x = x_i + radius*np.cos(t)
            ring_y = y_i + radius*np.sin(t)
            ring_z = np.full_like(t, z_i)
            ax.plot(ring_x, ring_y, ring_z, color=colors[i], alpha=0.25, linewidth=1.5)
        else:
            x_i, y_i = positions[i][:2]
            ax.plot(x_i, y_i, 'o', color=colors[i], label=f"{label_prefix}{i}")
            circle = plt.Circle((x_i, y_i), radius, color=colors[i], alpha=0.2)
            ax.add_patch(circle)

        for j in range(i + 1, n):
            if is3d:
                x_j, y_j, z_j = positions[j][:3]
                if np.linalg.norm(positions[j][:3] - positions[i][:3]) <= max_dist:
                    ax.plot([x_i, x_j], [y_i, y_j], [z_i, z_j], '-', color='gray', linewidth=1)
            else:
                x_j, y_j = positions[j][:2]
                if np.hypot(x_j - x_i, y_j - y_i) <= max_dist:
                    ax.plot([x_i, x_j], [y_i, y_j], '-', color='gray', linewidth=1)

def plot_trajectory(ax, traj, colors, is3d=False):
    """
    traj: (N, T, D)，D>=2；is3d=True 时取前三列 (x,y,z)，否则取前两列 (x,y)。
    """
    n = traj.shape[0]
    for i in range(n):
        if is3d:
            ax.plot(traj[i, :, 0], traj[i, :, 1], traj[i, :, 2], '-', color=colors[i], alpha=0.7)
        else:
            ax.plot(traj[i, :, 0], traj[i, :, 1], '-', color=colors[i], alpha=0.7)


def main():
    parser = argparse.ArgumentParser(description="Plot formation connectivity and trajectories")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--states", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--create_anim", action="store_true")
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

    HAS_RESULT = False
    if os.path.exists(states_file):
        st = load_json(states_file)
        HAS_RESULT = True

    # 初始和最终的机器人数组
    sf = np.array(cfg['tasks']['sf'], dtype=float)
    so = np.array(cfg['tasks']['so'], dtype=float)
    max_dist = cfg['cbf_params']['d_max']
    ts = cfg["mpc_params"]["h"]
    N = len(sf)

    # robots 数据结构：st['robots'][key]['states'] 是一个帧序列
    T = 0
    if HAS_RESULT:
        robots = st["robots"]
        keys = sorted(robots.keys(), key=int)  # 假设机器人 ID 是字符串数字
        traj = np.array([robots[k]["states"] for k in keys], dtype=float)  # (N, T, D)
        _, T, _ = traj.shape

    # 是否 3D：sf/so 或 traj 有第三列就启用 3D
    is3d = (sf.shape[1] >= 3) or (so.shape[1] >= 3)
    if HAS_RESULT:
        is3d = is3d or (traj.shape[2] >= 3)

    colors = generate_colors(N)

    # === 3. 创建子图 ===
    if HAS_RESULT:
        fig = plt.figure(figsize=(15, 5))
        axes = [
            fig.add_subplot(1, 3, 1, projection='3d' if is3d else None),
            fig.add_subplot(1, 3, 2, projection='3d' if is3d else None),
            fig.add_subplot(1, 3, 3, projection='3d' if is3d else None),
        ]
    else:
        fig = plt.figure(figsize=(10, 5))
        axes = [
            fig.add_subplot(1, 2, 1, projection='3d' if is3d else None),
            fig.add_subplot(1, 2, 2, projection='3d' if is3d else None),
        ]

    # --- 3.1 初始连通性（Static） ---
    axes[0].set_title('Initial Connectivity')
    plot_connectivity(axes[0], so, max_dist, colors, robot_radius, is3d=is3d)
    axes[0].grid(True)

    # --- 3.2 最终连通性（Static） ---
    axes[1].set_title('Final Connectivity')
    plot_connectivity(axes[1], sf, max_dist, colors, robot_radius, is3d=is3d)
    axes[1].grid(True)

    # --- 3.3 轨迹（仅当有结果时） ---
    if HAS_RESULT:
        axes[2].set_title("Trajectories & Connectivity Over Time")
        plot_trajectory(axes[2], traj, colors, is3d=is3d)
        axes[2].grid(True)

    # 轴标签与等比例
    for ax in axes:
        if is3d:
            ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
            try:
                ax.set_box_aspect([1, 1, 1])  # matplotlib>=3.3
            except Exception:
                pass

    # 如果配置里有物理边界，设置显示范围
    if "physical_limits" in cfg and "p_min" in cfg["physical_limits"] and "p_max" in cfg["physical_limits"]:
        p_min = np.array(cfg["physical_limits"]["p_min"], dtype=float)
        p_max = np.array(cfg["physical_limits"]["p_max"], dtype=float)
        for ax in axes:
            if is3d:
                ax.set_xlim(p_min[0], p_max[0])
                ax.set_ylim(p_min[1], p_max[1])
                if len(p_min) >= 3 and len(p_max) >= 3:
                    ax.set_zlim(p_min[2], p_max[2])
            else:
                ax.set_xlim(p_min[0], p_max[0])
                ax.set_ylim(p_min[1], p_max[1])

    # 保存静态图
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

            # 取出当前帧位置
            curr_positions = traj[:, frame_idx, :3] if is3d else traj[:, frame_idx, :2]

            # 1) 当前位置散点
            if is3d:
                scat = ax_anim.scatter(
                    curr_positions[:, 0],
                    curr_positions[:, 1],
                    curr_positions[:, 2],
                    c=colors,
                    s=50,
                    edgecolors='k',
                    zorder=3
                )
            else:
                scat = ax_anim.scatter(
                    curr_positions[:, 0],
                    curr_positions[:, 1],
                    c=colors,
                    s=50,
                    edgecolors='k',
                    zorder=3
                )
            artists.append(scat)

            # 2) 半径圈/环
            for i in range(N):
                if is3d:
                    x_i, y_i, z_i = curr_positions[i]
                    t = np.linspace(0, 2*np.pi, 80)
                    ring_x = x_i + robot_radius*np.cos(t)
                    ring_y = y_i + robot_radius*np.sin(t)
                    ring_z = np.full_like(t, z_i)
                    ring, = ax_anim.plot(ring_x, ring_y, ring_z, color=colors[i], alpha=0.25, zorder=1)
                    artists.append(ring)
                else:
                    x_i, y_i = curr_positions[i]
                    circle = plt.Circle((x_i, y_i), robot_radius, color=colors[i], alpha=0.2, zorder=1)
                    ax_anim.add_patch(circle)
                    artists.append(circle)

            # 3) 连通性连线
            for i in range(N):
                for j in range(i + 1, N):
                    if is3d:
                        if np.linalg.norm(curr_positions[j] - curr_positions[i]) <= max_dist:
                            ln, = ax_anim.plot(
                                [curr_positions[i, 0], curr_positions[j, 0]],
                                [curr_positions[i, 1], curr_positions[j, 1]],
                                [curr_positions[i, 2], curr_positions[j, 2]],
                                '-', color='gray', lw=1, zorder=2
                            )
                            artists.append(ln)
                    else:
                        dx, dy = curr_positions[j] - curr_positions[i]
                        if np.hypot(dx, dy) <= max_dist:
                            ln, = ax_anim.plot(
                                [curr_positions[i, 0], curr_positions[j, 0]],
                                [curr_positions[i, 1], curr_positions[j, 1]],
                                '-', color='gray', lw=1, zorder=2
                            )
                            artists.append(ln)

            return artists

        # Downsample frames if T is too large for smooth animation
        MAX_FRAMES = 200
        if T > MAX_FRAMES:
            downsample_factor = max(1, T // MAX_FRAMES)
            frame_indices = np.arange(0, T, downsample_factor)
            actual_frames = len(frame_indices)
            actual_interval = 1000 * ts * downsample_factor
            print(f"Downsampling animation: {T} frames -> {actual_frames} frames (factor: {downsample_factor})")

            def update_frame_downsampled(frame_idx):
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
