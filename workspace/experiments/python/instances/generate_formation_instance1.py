import json
import argparse
import os

def extract_max_dist(cpp_path: str) -> float:
    """直接返回固定的 max_dist=2.0，也可以根据cpp_path解析得到实际值"""
    return 2.0

def main(cpp_path: str, base_config_path: str, out_path: str):
    # 1. 读出 max_dist
    max_dist = extract_max_dist(cpp_path)

    # 2. 从 base_config.json 里读取配置
    if not os.path.isfile(base_config_path):
        raise FileNotFoundError(f"找不到 base_config 文件：{base_config_path}")
    with open(base_config_path, 'r', encoding='utf-8') as f:
        cfg = json.load(f)

    # 3. 用 extract_max_dist 返回值覆盖 connectivity_params.max_distance
    if "connectivity_params" not in cfg:
        cfg["connectivity_params"] = {}
    cfg["connectivity_params"]["max_distance"] = max_dist

    # 4. 确保输出目录存在
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)

    # 5. 写出最终 JSON
    with open(out_path, 'w', encoding='utf-8') as f_out:
        json.dump(cfg, f_out, indent=4)
    print(f"✅ 已生成配置：{out_path}，其中 max_distance = {max_dist}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="生成配置文件（从 base_config.json 读取并覆盖 max_distance）"
    )
    parser.add_argument(
        '-c', '--cpp',
        type=str,
        default='/usr/src/mpc-cbf/workspace/lib/cbf/examples/connectivity/CBFFormationControl_example.cpp',
        help='CBFFormationControl_example.cpp 路径，用于 extract_max_dist'
    )
    parser.add_argument(
        '-b', '--base',
        type=str,
        default='/usr/src/mpc-cbf/workspace/experiments/config/base_config.json',
        help='基础配置文件 base_config.json 的路径'
    )
    parser.add_argument(
        '-o', '--out',
        type=str,
        default='/usr/src/mpc-cbf/workspace/experiments/config/circle/circle3_config.json',
        help='输出的最终 JSON 文件路径'
    )
    args = parser.parse_args()
    main(args.cpp, args.base, args.out)
