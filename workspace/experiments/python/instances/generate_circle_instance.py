import argparse
import json
import os

import numpy as np
from utils import generate_points_on_circle


def main(out_path: str, num_robots: int, radius: float):
    
    start_x, start_y = generate_points_on_circle(num_robots, radius, 0)
    goal_x, goal_y = generate_points_on_circle(num_robots, radius, np.pi)


    output = {
        "tasks": {
            "so": [[start_x[i], start_y[i], 0] for i in range(num_robots)],
            "sf": [[goal_x[i], goal_y[i], 0] for i in range(num_robots)],
        }
    }
    
    # write to output_path
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    
    with open(out_path, 'w', encoding='utf-8') as f_out:
        json.dump(output, f_out, indent=4)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate circle instance configuration")
    parser.add_argument(
        '-o', '--out',
        type=str,
        help='Output JSON file path'
    )
    parser.add_argument(
        '-n', '--num_robots',
        type=int,
        required=True,
        help='Number of robots in the circle formation'
    )
    parser.add_argument(
        '-r', '--radius',
        type=float,
        default=2.0,
        help='Radius of the circle'
    )
    
    args = parser.parse_args()
    out = args.out or f'/usr/src/mpc-cbf/workspace/experiments/config/examples/{args.num_robots}r/circle.json'
    main(out, args.num_robots, args.radius)