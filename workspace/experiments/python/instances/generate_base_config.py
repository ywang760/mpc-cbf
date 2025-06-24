import json

import os

# 硬编码默认输出路径
output_path = '../../config/base_config.json'

# 确保输出目录存在
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# 基本配置字典
base_config = {
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
    "connectivity_params": {
        "max_distance": 0.0
    },
    "tasks": {
        "so": [
            [0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 3.0, 0.0]
        ],
        "sf": [
            [0.0, 1.0, 0.0],
            [1.5, 2.0, 0.0],
            [0.0, 3.0, 0.0]
        ]
    }
}

# 将配置写入 JSON 文件
with open(output_path, 'w') as f:
    json.dump(base_config, f, indent=4)

print(f"Generated base_config.json at: {output_path}")
