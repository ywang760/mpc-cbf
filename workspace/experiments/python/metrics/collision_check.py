import json

import numpy as np


def load_states(json_filename):
    f = open(json_filename)
    data = json.load(f)
    return data

def rectangles_collide(xA, yA, widthA, heightA, xB, yB, widthB, heightB):
    """
        here xA, yA is the left top corner of recA, widthA, heightA is the width and height of recA
    """
    return (
            xA < xB + widthB and
            xA + widthA > xB and
            yA < yB + heightB and
            yA + heightA > yB
    )

def collision_check(x1, y1, x2, y2, collision_shape):
    """x1, y1 is the center of pos, collision shape is the half of the collision width and height"""
    collision_x = collision_shape[0]
    collision_y = collision_shape[1]
    xA = x1 - collision_x/2
    yA = y1 - collision_y/2
    xB = x2 - collision_x/2
    yB = y2 - collision_y/2
    widthA = 2*collision_x
    heightA = 2*collision_y
    widthB = 2*collision_x
    heightB = 2*collision_y
    return rectangles_collide(xA, yA, widthA, heightA, xB, yB, widthB, heightB)

def reach_goal_area(pos, goal, radius=1):
    dist = np.linalg.norm(pos - goal)
    return dist <= radius

def instance_success(traj, goals, radius, collision_shape):
    """traj: [n_robot, ts, dim]"""
    n_robot = traj.shape[0]
    ts = traj.shape[1]
    reach_goals = [False for i in range(n_robot)]

    for i in range(n_robot):
        pos = traj[i, -1, :2]
        goal = goals[i][:2]
        if not reach_goal_area(pos, goal, radius):
            print("Cannot reach goal area...")
            # return False, float('inf')

    for t in range(ts):
        if all(reach_goals):
            return True, max(0, t-1)
        for i in range(n_robot):
            pos_1 = traj[i, t, :3]  # [3, ]
            goal_1 = goals[i]
            if reach_goal_area(pos_1[:2], goal_1[:2], radius):
                reach_goals[i] = True
            for j in range(i+1, n_robot):
                pos_2 = traj[j, t, :3]  # [3, ]
                if collision_check(pos_1[0], pos_1[1], pos_2[0], pos_2[1], collision_shape):
                    d = np.linalg.norm(pos_1[:2] - pos_2[:2])
                    print("collision happens at timestep: ", t)
                    print(
                        f"robot {i} at ({pos_1[0]}, {pos_1[1]}) collides with robot {j} at ({pos_2[0]}, {pos_2[1]}), distance: {d:.2f}"
                    )
                    return False, float('inf')

    return True, ts

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Check collision and success of robot trajectories.")
    parser.add_argument("--config", type=str, required=True, help="Path to the configuration JSON file.")
    parser.add_argument("--states", type=str, required=True, help="Path to the states JSON file.")
    args = parser.parse_args()
    config_json = load_states(args.config)
    states_json = load_states(args.states)
    num_robots = len(states_json["robots"])
    traj = np.array([states_json["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
    collision_shape = config_json["robot_params"]["collision_shape"]["aligned_box"][:2]
    goals = np.array(config_json["tasks"]["sf"])
    success = instance_success(traj, goals, radius=1, collision_shape=collision_shape)
