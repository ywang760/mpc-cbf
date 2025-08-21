import argparse
import json
import numpy as np
import os
import colorsys
import matplotlib
import matplotlib.animation as animation
matplotlib.use("Agg")  # 纯离屏渲染
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle

# ======================
# Utils
# ======================
def generate_rgb_colors(num_colors):
    output = []
    num_colors += 1  # avoid hsv=0
    for index in range(1, num_colors):
        output.append(colorsys.hsv_to_rgb(index / num_colors, 0.75, 0.75))
    return np.asarray(output)

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)

def _build_connectivity_segments(pts_xy, d_connect, color_mode="const"):
    """
    根据距离阈值 d_connect 生成连线段。
    color_mode: "const" | "dist_alpha" | "dist_gray"
    """
    N = pts_xy.shape[0]
    segments, colors = [], []
    for i in range(N):
        for j in range(i + 1, N):
            dx, dy = pts_xy[j] - pts_xy[i]
            dist = np.hypot(dx, dy)
            if d_connect is not None and dist <= d_connect:
                segments.append([pts_xy[i], pts_xy[j]])
                if color_mode == "dist_alpha":
                    a = max(0.2, 1.0 - dist / max(d_connect, 1e-9))
                    colors.append((0.3, 0.3, 0.3, a))
                elif color_mode == "dist_gray":
                    g = 0.2 + 0.8 * (dist / max(d_connect, 1e-9))
                    g = np.clip(g, 0.2, 1.0)
                    colors.append((g, g, g, 0.9))
    if color_mode == "const":
        return segments, None
    return segments, (np.array(colors) if len(colors) else None)


# ======================
# 三联图（单张横向三幅）
# ======================
def _draw_positions_panel(
    ax,
    positions,            # (N,2)
    colors,
    bbox_half=None,       # 备选：没半径时使用
    d_connect=None,
    title="",
    show_connectivity=True,
    conn_color_mode="dist_alpha",
    xlim=None, ylim=None,
    disk_radius=None      # 优先：半径>0 画圆圈
):
    N = positions.shape[0]
    ax.set_aspect("equal")
    ax.set_title(title)
    # ax.axis("off")   # 删掉这句，保留坐标轴/刻度
    ax.grid(True, alpha=0.25)  # 可选：更易比对


    # 统一坐标轴范围
    if xlim is not None and ylim is not None:
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    else:
        pad = 3.0
        ax.set_xlim(positions[:,0].min()-pad, positions[:,0].max()+pad)
        ax.set_ylim(positions[:,1].min()-pad, positions[:,1].max()+pad)

    # 位置点
    for i in range(N):
        ax.plot(positions[i,0], positions[i,1],
        marker='o', ms=6, ls='None', c=colors[i], alpha=0.95, zorder=2)


    # 形状：优先圆圈，其次矩形
    if disk_radius is not None and disk_radius > 0:
        for i in range(N):
            circ = Circle(
                (positions[i,0], positions[i,1]),
                radius=disk_radius,
                ec=colors[i],           # 边线颜色：机器人专属色
                fc=colors[i],           # 填充颜色：机器人专属色
                lw=1.2,
                linestyle="--",         # 虚线边框
                alpha=0.25,             # 填充透明度（看得见连接线/轨迹）
                zorder=1,               # 圆圈放底层
            )
            ax.add_patch(circ)

    elif bbox_half is not None:
        for i in range(N):
            rect = plt.Rectangle((positions[i,0]-bbox_half[0], positions[i,1]-bbox_half[1]),
                                 2*bbox_half[0], 2*bbox_half[1],
                                 lw=1.0, ec='blue', fc='none', alpha=0.85)
            ax.add_patch(rect)

    # 连通性
    if show_connectivity and d_connect is not None:
        segs, seg_colors = _build_connectivity_segments(positions, d_connect, color_mode=conn_color_mode)
        if segs:
            lc = LineCollection(segs,
                    colors=seg_colors if seg_colors is not None else 'gray',
                    lw=1.2, alpha=0.9, zorder=3)

            ax.add_collection(lc)

def save_triptych(config, traj, colors, out_path, d_connect,
                  conn_color_mode="dist_alpha"):
    """
    生成单张横向三联图：
      左：初始 so
      中：目标 sf
      右：实际最终帧
    三幅共享同一坐标轴范围（若 physical_limits 给出）。
    机器人“范围”优先使用 config.collision_shape.radius 的圆圈，否则回退到 aligned_box 矩形。
    """
    so = np.array(config["tasks"]["so"], dtype=float)[:, :2]
    sf = np.array(config["tasks"]["sf"], dtype=float)[:, :2]
    final = traj[:, -1, :2]

    # 半径（优先）
    disk_radius = None
    try:
        disk_radius = float(config["robot_params"]["collision_shape"]["radius"])
    except Exception:
        disk_radius = None

    # 矩形半尺寸（备选）
    bbox_half = None
    try:
        bbox_half = np.array(config["robot_params"]["collision_shape"]["aligned_box"][:2], dtype=float)
    except Exception:
        bbox_half = None

    # 三联图统一轴范围：若 physical_limits 存在则用它，否则用三帧联合的 min/max
    xlim = ylim = None
    if "physical_limits" in config and "p_min" in config["physical_limits"] and "p_max" in config["physical_limits"]:
        p_min = np.array(config["physical_limits"]["p_min"], dtype=float)
        p_max = np.array(config["physical_limits"]["p_max"], dtype=float)
        xlim = (p_min[0], p_max[0]); ylim = (p_min[1], p_max[1])
    else:
        all_pts = np.vstack([so, sf, final])
        pad = 3.0
        xlim = (all_pts[:,0].min()-pad, all_pts[:,0].max()+pad)
        ylim = (all_pts[:,1].min()-pad, all_pts[:,1].max()+pad)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)

    _draw_positions_panel(axs[0], so,    colors, bbox_half, d_connect,
                          title="Initial (so)", show_connectivity=True,
                          conn_color_mode=conn_color_mode,
                          xlim=xlim, ylim=ylim, disk_radius=disk_radius)
    _draw_positions_panel(axs[1], sf,    colors, bbox_half, d_connect,
                          title="Desired (sf)", show_connectivity=True,
                          conn_color_mode=conn_color_mode,
                          xlim=xlim, ylim=ylim, disk_radius=disk_radius)
    _draw_positions_panel(axs[2], final, colors, bbox_half, d_connect,
                          title="Actual Final", show_connectivity=True,
                          conn_color_mode=conn_color_mode,
                          xlim=xlim, ylim=ylim, disk_radius=disk_radius)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Triptych saved: {out_path}")


# ======================
# 可选：动画（同样使用 radius 画圆圈，缺失时回退矩形）
# ======================
def animation2D_XYYaw(
    traj,
    dt,
    Ts,
    bbox_half,
    pred_curve=None,
    Trailing=True,
    PredCurve=True,
    colors="r",
    save_name="./test.mp4",
    d_connect=None,
    show_connectivity=True,
    conn_color_mode="const",
    disk_radius=None,   # ⬅️ 新增：圆圈半径（来自 config["robot_params"]["collision_shape"]["radius"]）
):
    """
    动画：彩色填充虚线圆圈 + 位置点 + 速度向量 + 尾迹 + 连通性连线（无 λ₂）。
    """
    n_agent, total_frame, _ = traj.shape
    trailing_buf_size = 1000000
    dir_len = 0.1

    x = traj[:, :, 0]
    y = traj[:, :, 1]
    x_v = traj[:, :, 3]
    y_v = traj[:, :, 4]

    x_min = np.min(x) - 3
    x_max = np.max(x) + 3
    y_min = np.min(y) - 3
    y_max = np.max(y) + 3
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.25)  # 保留坐标轴和网格，便于对比

    # ====== 初始帧元素 ======

    # 圆圈（彩色填充 + 虚线边）
    disks = []
    if disk_radius is not None and disk_radius > 0:
        for i in range(n_agent):
            circ = Circle(
                (x[i, 0], y[i, 0]),
                radius=disk_radius,
                ec=colors[i],
                fc=colors[i],
                lw=1.2,
                linestyle="--",
                alpha=0.25,   # 填充透明
                zorder=1      # 圆圈在底层
            )
            ax.add_patch(circ)
            disks.append(circ)

    # 点（在圆圈之上）
    p = [ax.plot([x[i, 0]], [y[i, 0]],
                 c=colors[i], marker='o', ms=6, ls='None', alpha=0.95, zorder=2)
         for i in range(n_agent)]

    # 速度线（点之上/连通性之下）
    vel_line = []
    for i in range(n_agent):
        (ln,) = ax.plot([x[i, 0], x[i, 0] + dir_len * x_v[i, 0]],
                        [y[i, 0], y[i, 0] + dir_len * y_v[i, 0]],
                        c='dodgerblue', alpha=0.7, zorder=2)
        vel_line.append(ln)

    # 尾迹（放在线下/点上）
    trail = []
    if Trailing:
        for i in range(n_agent):
            (ln,) = ax.plot([x[i, 0]], [y[i, 0]], c=colors[i], lw=3, alpha=0.4, zorder=2)
            trail.append(ln)

    # 预测曲线（若有）
    preds = None
    # 预测曲线（若有）
    preds = None
    if PredCurve and pred_curve is not None:
        impc_iters = pred_curve.shape[-3]
        preds = [[None for _ in range(impc_iters)] for _ in range(n_agent)]
        for i in range(n_agent):
            for it in range(impc_iters):
                (ln,) = ax.plot([], [], 
                                c=colors[i],       # ⬅️ 改成跟该 agent 一致的颜色
                                alpha=0.8,
                                linestyle="-",     # 预测曲线保持实线
                                zorder=2)
                preds[i][it] = ln

    # 连通性线（放在最上层）
    conn_lc = None
    if show_connectivity and d_connect is not None:
        conn_lc = LineCollection([], colors='gray', lw=1.2, alpha=0.9, zorder=3)
        ax.add_collection(conn_lc)

    # ====== 帧更新函数 ======
    def animate(t):
        updated = []

        # 圆圈中心更新
        if disks:
            for i in range(n_agent):
                disks[i].center = (x[i, t], y[i, t])
            updated += disks

        # 点 & 速度
        for i in range(n_agent):
            p[i][0].set_data([x[i, t]], [y[i, t]])
            vel_line[i].set_data([x[i, t], x[i, t] + dir_len * x_v[i, t]],
                                 [y[i, t], y[i, t] + dir_len * y_v[i, t]])
        updated += [p[i][0] for i in range(n_agent)] + vel_line

        # 尾迹
        if Trailing:
            tail = t if t < trailing_buf_size else trailing_buf_size
            for i in range(n_agent):
                trail[i].set_data(x[i, t - tail:t + 1], y[i, t - tail:t + 1])
            updated += trail

        # 预测曲线
        if PredCurve and pred_curve is not None:
            impc_iters = pred_curve.shape[-3]
            pred_idx = int(t // max(int(dt / Ts), 1))
            for i in range(n_agent):
                for it in range(impc_iters):
                    curve = pred_curve[:, :, it, :, :]
                    preds[i][it].set_data(curve[i, pred_idx, :, 0], curve[i, pred_idx, :, 1])
            updated += [preds[i][it] for i in range(n_agent) for it in range(impc_iters)]

        # 连通性
        if conn_lc is not None:
            pts = np.column_stack((x[:, t], y[:, t]))
            segs, seg_colors = _build_connectivity_segments(pts, d_connect, color_mode=conn_color_mode)
            conn_lc.set_segments(segs)
            conn_lc.set_color(seg_colors if seg_colors is not None else 'gray')
            updated.append(conn_lc)

        return tuple(updated)

    # 生成并保存
    anim = animation.FuncAnimation(fig, animate, frames=total_frame, interval=Ts * 1e3, blit=True)
    os.makedirs(os.path.dirname(save_name) or ".", exist_ok=True)
    try:
        anim.save(save_name, writer=animation.FFMpegWriter(fps=1 / Ts))
        print(f"[viz] Animation saved: {save_name}")
    except Exception as e:
        print(f"[viz] Failed to create/save animation: {e}")
    return anim


# ======================
# Main
# ======================
def main():
    print(">>> entered main")
    parser = argparse.ArgumentParser(description="Triptych visualization (radius circles) + optional animation")
    parser.add_argument("--config", type=str, required=True, help="path to config json file")
    parser.add_argument("--states", type=str, required=True, help="path to simulation state json file")
    parser.add_argument("--output_dir", type=str, required=True, help="directory to save figures/animations")
    parser.add_argument("--create_anim", action="store_true", help="create animation (mp4)")
    parser.add_argument("--anim_format", type=str, choices=["mp4", "gif"], default="gif",
                        help="Kept for compatibility; mp4 will be used.")
    args = parser.parse_args()

    cfg = load_json(args.config)
    if not os.path.exists(args.states):
        raise FileNotFoundError(f"States file not found: {args.states}")
    st = load_json(args.states)
    disk_radius = None
    try:
        disk_radius = float(cfg["robot_params"]["collision_shape"]["radius"])
    except Exception:
        disk_radius = None

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

    # 轨迹 (N, T, D)
    keys = sorted(st["robots"].keys(), key=int)
    traj = np.array([st["robots"][k]["states"] for k in keys], dtype=float)

    # 时间参数
    dt = st["dt"]
    Ts = st["Ts"]

    # 备选的矩形半尺寸（动画中若无半径使用）
    bbox_half = None
    try:
        bbox_half = np.array(cfg["robot_params"]["collision_shape"]["aligned_box"][:2], dtype=float)
    except Exception:
        bbox_half = None

    # 连通性阈值
    d_connect = cfg.get("cbf_params", {}).get("d_max", None)

    # === 三联图（圆圈范围优先） ===
    save_triptych(cfg, traj, colors, triptych_path, d_connect,
                  conn_color_mode="dist_alpha")

    # === 动画（可选，同样优先圆圈） ===
    if args.create_anim:
        pred_curve = None
        if "pred_curve" in st["robots"]["0"]:
            pred_curve = np.array([st["robots"][str(_)]["pred_curve"] for _ in range(num_robots)])

        # 读取半径（动画也用）
        disk_radius = None
        try:
            disk_radius = float(cfg["robot_params"]["collision_shape"]["radius"])
        except Exception:
            disk_radius = None

        animation2D_XYYaw(
            traj=traj,
            dt=dt,
            Ts=Ts,
            bbox_half=bbox_half,
            pred_curve=pred_curve,
            Trailing=True,
            PredCurve=(pred_curve is not None),
            colors=colors,
            save_name=anim_path,
            d_connect=d_connect,
            show_connectivity=(d_connect is not None),
            conn_color_mode="dist_alpha",
            disk_radius=disk_radius,  # ⬅️ 新增
        )
    else:
        print("[viz] Animation creation skipped.")

if __name__ == "__main__":
    main()
