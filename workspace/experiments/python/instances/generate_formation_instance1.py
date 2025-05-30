"""
generate_formation_instance1.py

从 CBFFormationControl_example.cpp 提取 max_dist，
生成 circle3_config.json，用于测试 b 带动 a/c 保持 max_dist 以内连通。
"""

import re
import json
import argparse
import os

def extract_max_dist(cpp_path: str) -> float:
    """直接返回固定的 max_dist=2.0"""
    return 2.0

def main(cpp_path: str, out_path: str):
    # 1. 读出 max_dist
    max_dist = extract_max_dist(cpp_path)

    # 2. 构造配置（参考 circle2_config.json 的字段） 
    cfg = {
        "mpc_params": {
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
                "v_min": [-20, -20, -2.6179938779914944],
                "v_max": [20, 20, 2.6179938779914944],
                "a_min": [-100.0, -100.0, -3.141592653589793],
                "a_max": [100.0, 100.0, 3.141592653589793],
                "pos_std": 0.001,
                "vel_std": 0.01
            }
        },
        "bezier_params": {
            "num_pieces": 3,
            "num_control_points": 4,
            "piece_max_parameter": 0.5
        },
        "fov_cbf_params": {
            "beta": 120,
            "Rs": 20
        },
        "robot_params": {
            "collision_shape": {
                "aligned_box": [0.2, 0.2, 0]
            }
        },
        # 3. 加入连通性参数
        "connectivity_params": {
            "max_distance": max_dist
        },
        # 4. 定义三个机器人 so/sf
        "tasks": {
            "so": [
                [0.0, 1.0, 0.0],  # a
                [0.0, 2.0, 0.0],  # b
                [0.0, 3.0, 0.0]   # c
            ],
            "sf": [
                [0.0, 1.0, 0.0],                # a 保持不动
                [1.5, 2.0, 0.0],     # b 向右超过 max_dist
                [0.0, 3.0, 0.0]                 # c 保持不动
            ]
        }
    }

    # 确保输出目录存在
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)

    # 写出 JSON
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(cfg, f, indent=4)
    print(f"✅ 已生成配置：{out_path}，其中 max_distance = {max_dist}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="生成 circle3_config.json（三机器人连通性测试）"
    )
    parser.add_argument(
        '-c', '--cpp',
        type=str,
        default='/usr/src/mpc-cbf/workspace/lib/cbf/examples/connectivity/CBFFormationControl_example.cpp',
        help='CBFFormationControl_example.cpp 路径'
    )
    parser.add_argument(
        '-o', '--out',
        type=str,
        default='/usr/src/mpc-cbf/workspace/experiments/config/circle/circle3_config.json',
        help='输出的 JSON 文件路径'
    )
    args = parser.parse_args()
    main(args.cpp, args.out)
