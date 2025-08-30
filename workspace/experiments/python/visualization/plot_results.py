import argparse
import json
import numpy as np
import os
import colorsys
import matplotlib
matplotlib.use("Agg")  # 纯离屏渲染（服务器/无显示环境）
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import art3d
import matplotlib.animation as animation


# ======================
# Utils
# ======================
def generate_rgb_colors(num_colors: int):
    """生成区分度较高的 RGB 颜色组。"""
    output = []
    num_colors += 1  # avoid hsv=0
    for index in range(1, num_colors):
        output.append(colorsys.hsv_to_rgb(index / num_colors, 0.75, 0.75))
    return np.asarray(output)

def _darker(color, factor=0.7):
    """把颜色加深：factor 越小越深 (0~1)。"""
    return tuple(np.clip(np.array(color[:3]) * factor, 0, 1))

def load_json(path: str):
    with open(path, "r") as f:
        return json.load(f)

# 可选：等比例 3D 轴（未默认启用）
def _set_axes_equal(ax):
    """让 3D 坐标轴比例一致（避免被压扁）。"""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0]); x_mid = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0]); y_mid = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0]); z_mid = np.mean(z_limits)
    r = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_mid - r, x_mid + r])
    ax.set_ylim3d([y_mid - r, y_mid + r])
    ax.set_zlim3d([z_mid - r, z_mid + r])

def _compute_bounds(xyz, cfg, pad=1.0, include_so_sf=True, disk_radius=None):
    """返回 x/y/z 的 [min, max]，包含整个轨迹，可选叠加 so/sf 与圆盘半径。"""
    x_all = xyz[:, :, 0]; y_all = xyz[:, :, 1]; z_all = xyz[:, :, 2]
    xmin, xmax = np.min(x_all), np.max(x_all)
    ymin, ymax = np.min(y_all), np.max(y_all)
    zmin, zmax = np.min(z_all), np.max(z_all)

    if include_so_sf:
        try:
            so = np.array(cfg["tasks"]["so"], dtype=float)
            sf = np.array(cfg["tasks"]["sf"], dtype=float)
            xmin = min(xmin, so[:,0].min(), sf[:,0].min())
            xmax = max(xmax, so[:,0].max(), sf[:,0].max())
            ymin = min(ymin, so[:,1].min(), sf[:,1].min())
            ymax = max(ymax, so[:,1].max(), sf[:,1].max())
            zmin = min(zmin, so[:,2].min(), sf[:,2].min())
            zmax = max(zmax, so[:,2].max(), sf[:,2].max())
        except Exception:
            pass

    # 圆盘半径留出边距
    rpad = float(disk_radius) if (disk_radius is not None) else 0.0
    xmin -= (pad + rpad); xmax += (pad + rpad)
    ymin -= (pad + rpad); ymax += (pad + rpad)
    if zmin == zmax:
        zmin -= 1.0; zmax += 1.0
    zmin -= pad; zmax += pad
    return [xmin, xmax], [ymin, ymax], [zmin, zmax]


# ======================
# 保存“任意帧”的静态图片
# ======================
def save_snapshot_at_frame(
    cfg,
    traj,               # (N, T, D)
    colors,             # (N, 3) or (N, 4)
    save_path,          # 输出路径 *.png / *.jpg / *.pdf
    d_connect=None,
    conn_color_mode="dist_alpha",   # "const" | "dist_alpha" | "dist_gray"
    disk_radius=None,
    disks_at_agent_z=True,          # True: 圆盘画在各自 z；False: 画在 z=0
    frame_idx=-1,                   # -1: 终点；0: 初始；也可传任意 0..T-1
    show_trail=True,
    trail_alpha=0.6,
    trail_lw=2.0,
    axis_mode="fit_all",            # 坐标轴模式
    show_goals=True,                # 是否绘制目标点 X
):
    n_agent, total_frame, D = traj.shape
    xyz = traj[:, :, :3] if D >= 3 else np.pad(traj[:, :, :2], ((0,0),(0,0),(0,1)))
    frame_idx = int(frame_idx) % total_frame
    cur_xyz = xyz[:, frame_idx, :]     # (N, 3)

    # —— 坐标轴范围 ——
    if axis_mode == "fit_all":
        xlim, ylim, zlim = _compute_bounds(xyz, cfg, pad=1.0, include_so_sf=True, disk_radius=disk_radius)
        fig = plt.figure(figsize=(8,8)); ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlim)
    elif axis_mode == "auto":
        fig = plt.figure(figsize=(8,8)); ax = fig.add_subplot(111, projection="3d")
        ax.set_autoscale_on(True)
    else:  # fit_current
        fig = plt.figure(figsize=(8,8)); ax = fig.add_subplot(111, projection="3d")
        # 先给全局范围，随后再按当前帧对焦一次（见文末）
        xlim, ylim, zlim = _compute_bounds(xyz, cfg, pad=1.0, include_so_sf=True, disk_radius=disk_radius)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlim)

    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.grid(True, alpha=0.25)

    # 1) 已运动轨迹
    if show_trail and frame_idx > 0:
        for i in range(n_agent):
            xs = xyz[i, :frame_idx+1, 0]
            ys = xyz[i, :frame_idx+1, 1]
            zs = xyz[i, :frame_idx+1, 2]
            ax.plot(xs, ys, zs, c=colors[i], lw=trail_lw, alpha=trail_alpha)

    # 2) 圆盘（可选）
    if (disk_radius is not None) and (disk_radius > 0):
        for i in range(n_agent):
            edge = tuple(colors[i][:3]) if len(colors[i]) >= 3 else (0.3, 0.3, 0.3)
            face = edge + (0.55,)
            circ = Circle((cur_xyz[i, 0], cur_xyz[i, 1]), radius=float(disk_radius))
            circ.set_edgecolor(edge)
            circ.set_facecolor(face)
            circ.set_linewidth(1.8)
            circ.set_linestyle("--")
            ax.add_patch(circ)
            z_plane = float(cur_xyz[i, 2]) if disks_at_agent_z else 0.0
            art3d.pathpatch_2d_to_3d(circ, z=z_plane, zdir="z")

    # 3) 当前帧散点（机器人当前位姿）
    ax.scatter(cur_xyz[:, 0], cur_xyz[:, 1], cur_xyz[:, 2], s=40, c=colors, depthshade=True)

    # 4) 目标点（sf）标注（深一号颜色）
    if show_goals:
        try:
            sf = np.array(cfg["tasks"]["sf"], dtype=float)  # 形状 (N, 3)
            for i in range(min(n_agent, sf.shape[0])):
                dark_c = _darker(colors[i], factor=0.7)
                ax.scatter(
                    [sf[i, 0]], [sf[i, 1]], [sf[i, 2]],
                    s=180, marker="X", c=[dark_c],
                    edgecolors="k", linewidths=0.8, alpha=0.95, depthshade=False
                )
        except Exception:
            pass

    # 5) 连通边
    if (d_connect is not None):
        d_connect = float(d_connect)
        for i in range(n_agent):
            pi = cur_xyz[i]
            for j in range(i + 1, n_agent):
                pj = cur_xyz[j]
                dist = np.hypot(pj[0] - pi[0], pj[1] - pi[1])  # xy 距离；如需 3D：np.linalg.norm(pj-pi)
                if dist <= d_connect:
                    if conn_color_mode == "dist_alpha":
                        a = max(0.2, 1.0 - dist / max(d_connect, 1e-9))
                        color = (0.3, 0.3, 0.3, a)
                    elif conn_color_mode == "dist_gray":
                        g = 0.2 + 0.8 * (dist / max(d_connect, 1e-9))
                        g = np.clip(g, 0.2, 1.0)
                        color = (g, g, g, 0.9)
                    else:
                        color = "gray"
                    ax.plot([pi[0], pj[0]], [pi[1], pj[1]], [pi[2], pj[2]], c=color, lw=1.2)

    # 若选择 fit_current：按当前帧再对焦一次，确保取景跟随
    if axis_mode == "fit_current":
        pad = 0.5
        xmin = cur_xyz[:,0].min()-pad; xmax = cur_xyz[:,0].max()+pad
        ymin = cur_xyz[:,1].min()-pad; ymax = cur_xyz[:,1].max()+pad
        zmin = cur_xyz[:,2].min()-pad; zmax = cur_xyz[:,2].max()+pad
        if disk_radius:
            xmin -= disk_radius; xmax += disk_radius
            ymin -= disk_radius; ymax += disk_radius
        if zmin == zmax:
            zmin -= 1.0; zmax += 1.0
        ax.set_xlim([xmin, xmax]); ax.set_ylim([ymin, ymax]); ax.set_zlim([zmin, zmax])

    os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
    fig.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Snapshot saved: {save_path}")


# ======================
# 动画（方案 A：每帧重建圆盘）
# ======================
def animation3D_XYYaw(
    cfg,
    traj,
    dt,
    Ts,
    pred_curve=None,
    Trailing=True,
    PredCurve=True,
    colors="r",
    save_name="./test3d.mp4",
    d_connect=None,
    show_connectivity=True,
    conn_color_mode="const",
    disk_radius=None,
    disks_at_agent_z=True,   # True: 圆盘画在各自 z；False: 画在 z=0
    axis_mode="fit_all",     # 坐标轴模式
    show_goals=True,         # 是否绘制目标点 X
):
    n_agent, total_frame, D = traj.shape
    trailing_buf_size = 1000000
    dir_len = 0.1

    # 轨迹坐标
    x = traj[:, :, 0]
    y = traj[:, :, 1]
    z = traj[:, :, 2] if D > 2 else np.zeros_like(x)

    # 速度（若不存在则用 0）
    x_v = traj[:, :, 3] if D > 3 else np.zeros_like(x)
    y_v = traj[:, :, 4] if D > 4 else np.zeros_like(y)

    # —— 坐标轴范围 ——
    if axis_mode == "fit_all":
        xyz = traj[:, :, :3] if D >= 3 else np.pad(traj[:, :, :2], ((0,0),(0,0),(0,1)))
        xlim, ylim, zlim = _compute_bounds(xyz, cfg, pad=1.0, include_so_sf=True, disk_radius=disk_radius)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlim)
    elif axis_mode == "auto":
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_autoscale_on(True)
    else:  # fit_current
        xyz = traj[:, :, :3] if D >= 3 else np.pad(traj[:, :, :2], ((0,0),(0,0),(0,1)))
        xlim, ylim, zlim = _compute_bounds(xyz, cfg, pad=1.0, include_so_sf=True, disk_radius=disk_radius)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlim)

    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.grid(True, alpha=0.25)

    # 目标点（sf）—— 深一号颜色，动画中用稍小 s 以免遮挡
    if show_goals:
        try:
            sf = np.array(cfg["tasks"]["sf"], dtype=float)
            for i in range(min(n_agent, sf.shape[0])):
                dark_c = _darker(colors[i], factor=0.7)
                ax.scatter(
                    [sf[i, 0]], [sf[i, 1]], [sf[i, 2]],
                    s=90, marker="X", c=[dark_c],
                    edgecolors="k", linewidths=0.8, alpha=0.95, depthshade=False
                )
        except Exception:
            pass

    # 连接线对象（每对 agent 预创建一条）
    pairs = [(i, j) for i in range(n_agent) for j in range(i + 1, n_agent)]
    conn_lines = []
    if show_connectivity and d_connect is not None:
        for _ in pairs:
            (ln,) = ax.plot([], [], [], c='gray', lw=1.2, alpha=0.9)
            ln.set_visible(False)
            conn_lines.append(ln)

    # 点（当前位姿）
    points = []
    for i in range(n_agent):
        (ln,) = ax.plot([x[i, 0]], [y[i, 0]], [z[i, 0]], marker='o', ls='None',
                        c=colors[i], alpha=0.95)
        points.append(ln)

    # 速度向量
    vel_lines = []
    for i in range(n_agent):
        (ln,) = ax.plot([x[i, 0], x[i, 0] + dir_len * x_v[i, 0]],
                        [y[i, 0], y[i, 0] + dir_len * y_v[i, 0]],
                        [z[i, 0], z[i, 0]],
                        c='dodgerblue', alpha=0.7)
        vel_lines.append(ln)

    # 尾迹
    trails = []
    if Trailing:
        for i in range(n_agent):
            (ln,) = ax.plot([x[i, 0]], [y[i, 0]], [z[i, 0]], c=colors[i], lw=2, alpha=0.6)
            trails.append(ln)

    # —— 方案 A：每帧重建圆盘 ——
    live_disks = [None] * n_agent

    def animate(t):
        updated = []

        # 每帧重建圆盘（先移除上一帧的，再新建这一帧）
        if disk_radius is not None and disk_radius > 0:
            for i in range(n_agent):
                if live_disks[i] is not None:
                    try:
                        live_disks[i].remove()
                    except Exception:
                        pass
                    live_disks[i] = None

                edge = tuple(colors[i][:3]) if len(colors[i]) >= 3 else (0.3, 0.3, 0.3)
                face = edge + (0.55,)
                circ = Circle((x[i, t], y[i, t]), radius=float(disk_radius))
                circ.set_edgecolor(edge)
                circ.set_facecolor(face)
                circ.set_linewidth(1.8)
                circ.set_linestyle("--")
                ax.add_patch(circ)
                z_plane = float(z[i, t]) if disks_at_agent_z else 0.0
                art3d.pathpatch_2d_to_3d(circ, z=z_plane, zdir="z")
                live_disks[i] = circ
                updated.append(circ)

        # 点 & 速度线
        for i in range(n_agent):
            points[i].set_data_3d([x[i, t]], [y[i, t]], [z[i, t]])
            vel_lines[i].set_data_3d([x[i, t], x[i, t] + dir_len * x_v[i, t]],
                                     [y[i, t], y[i, t] + dir_len * y_v[i, t]],
                                     [z[i, t], z[i, t]])
        updated += points + vel_lines

        # 尾迹
        if Trailing:
            tail = t if t < trailing_buf_size else trailing_buf_size
            for i in range(n_agent):
                xs = x[i, max(0, t - tail):t + 1]
                ys = y[i, max(0, t - tail):t + 1]
                zs = z[i, max(0, t - tail):t + 1]
                trails[i].set_data_3d(xs, ys, zs)
            updated += trails

        # 连通性（复用对象，不再新建）
        if show_connectivity and d_connect is not None:
            d_connect_f = float(d_connect)
            for k, (i, j) in enumerate(pairs):
                dx = x[j, t] - x[i, t]
                dy = y[j, t] - y[i, t]
                dist = np.hypot(dx, dy)  # 如需 3D 距离：np.linalg.norm([dx, dy, z[j,t]-z[i,t]])
                ln = conn_lines[k]
                if dist <= d_connect_f:
                    ln.set_data_3d([x[i, t], x[j, t]],
                                   [y[i, t], y[j, t]],
                                   [z[i, t], z[j, t]])
                    if conn_color_mode == "dist_alpha":
                        a = max(0.2, 1.0 - dist / max(d_connect_f, 1e-9))
                        ln.set_alpha(a); ln.set_color((0.3, 0.3, 0.3))
                    elif conn_color_mode == "dist_gray":
                        g = 0.2 + 0.8 * (dist / max(d_connect_f, 1e-9))
                        g = np.clip(g, 0.2, 1.0)
                        ln.set_alpha(0.9); ln.set_color((g, g, g))
                    else:
                        ln.set_alpha(0.9); ln.set_color("gray")
                    ln.set_visible(True)
                else:
                    ln.set_visible(False)
            updated += conn_lines

        # 若选择 fit_current：按当前帧对焦
        if axis_mode == "fit_current":
            pad = 0.5
            xmin = x[:,t].min()-pad; xmax = x[:,t].max()+pad
            ymin = y[:,t].min()-pad; ymax = y[:,t].max()+pad
            zmin = z[:,t].min()-pad; zmax = z[:,t].max()+pad
            if disk_radius:
                xmin -= disk_radius; xmax += disk_radius
                ymin -= disk_radius; ymax += disk_radius
            if zmin == zmax:
                zmin -= 1.0; zmax += 1.0
            ax.set_xlim([xmin, xmax]); ax.set_ylim([ymin, ymax]); ax.set_zlim([zmin, zmax])

        return tuple(updated)

    # 帧序列 & 动画
    frames = np.arange(total_frame)
    anim = animation.FuncAnimation(fig, animate, frames=frames, interval=Ts * 1e3, blit=False)

    # 控制 fps，避免过高导致压制后看不清细节
    target_fps = min(30, max(5, int(1.0 / Ts)))
    writer = animation.FFMpegWriter(fps=target_fps)

    os.makedirs(os.path.dirname(save_name) or ".", exist_ok=True)
    try:
        anim.save(save_name, writer=writer)
        print(f"[viz] 3D animation saved: {save_name}")
    except Exception as e:
        print(f"[viz] Failed to create/save 3D animation: {e}")
    return anim


# ======================
# Main
# ======================
def main():
    print(">>> entered main")
    parser = argparse.ArgumentParser(description="Triptych visualization (true-3D) + optional animation")
    parser.add_argument("--config", type=str, required=True, help="path to config json file")
    parser.add_argument("--states", type=str, required=True, help="path to simulation state json file")
    parser.add_argument("--output_dir", type=str, required=True, help="directory to save figures/animations")
    parser.add_argument("--create_anim", action="store_true", help="create animation (mp4)")
    parser.add_argument("--anim_format", type=str, choices=["mp4", "gif"], default="gif",
                        help="Kept for compatibility; mp4 will be used.")
    parser.add_argument("--axis_mode", type=str, choices=["fit_all","fit_current","auto"], default="fit_all",
                        help="3D轴缩放模式：fit_all(全局固定)/fit_current(逐帧跟随)/auto(由Matplotlib自动)")
    parser.add_argument("--no_goals", action="store_true", help="不绘制终点 X 标记")  # 新增开关
    args = parser.parse_args()

    cfg = load_json(args.config)
    if not os.path.exists(args.states):
        raise FileNotFoundError(f"States file not found: {args.states}")
    st = load_json(args.states)

    # 输出路径
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    out_dir_expanded = os.path.join(out_dir, os.path.dirname(args.config).split("/")[-1])
    os.makedirs(out_dir_expanded, exist_ok=True)
    triptych_path = os.path.join(out_dir_expanded, os.path.basename(args.config).replace(".json", "_triptych.png"))
    anim_path     = os.path.join(out_dir_expanded, os.path.basename(args.config).replace(".json", ".mp4"))

    # 颜色
    num_robots = len(st["robots"])
    colors = generate_rgb_colors(num_robots)

    # 轨迹 (N, T, D) —— 按 robots 的 key 排序，保证与 so/sf 映射一致
    keys = sorted(st["robots"].keys(), key=int)
    traj = np.array([np.array(st["robots"][k]["states"], dtype=float) for k in keys])

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

    # 半径（从 config 读取）
    try:
        disk_radius = float(cfg["robot_params"]["collision_shape"]["radius"])
    except Exception:
        disk_radius = None

    # 连通性阈值
    d_connect = cfg.get("cbf_params", {}).get("d_max", None)

    # 终点显示与否
    show_goals = (not args.no_goals)

    # === 静态快照：多帧采样 ===
    total_frame = traj.shape[1]
    sample_frames = np.linspace(0, total_frame - 1, 6, dtype=int)  # 6个采样点

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

    # === 动画（可选） ===
    if args.create_anim:
        animation3D_XYYaw(
            cfg=cfg,
            traj=traj,
            dt=dt,
            Ts=Ts,
            pred_curve=None,
            Trailing=True,
            PredCurve=False,
            colors=colors,
            save_name=anim_path,
            d_connect=d_connect,
            show_connectivity=(d_connect is not None),
            conn_color_mode="dist_alpha",
            disk_radius=disk_radius,       # 半径来自 config
            disks_at_agent_z=True,         # 圆盘跟随各自 z 高度
            axis_mode=args.axis_mode,
            show_goals=show_goals,
        )
    else:
        print("[viz] Animation creation skipped.")

    print("All configurations processed successfully!")


if __name__ == "__main__":
    main()
