from Metrics import *

result_path = "../instances/results/log"
config_path = "../instances/"
date = "01012025"
instance_type = "circle"
min_r = 2
max_r = 8
num_robot = range(min_r, max_r+1)
num_exp = 3

success_rate = []  # [entry,], averaged by M, M is the sample
num_neighbor_in_fov = []  # [entry, M]

for num_r in num_robot:
    for exp_idx in range(3):
        state_json = result_path+date+"/"+instance_type+str(num_r)+"States_"+str(exp_idx)+".json"
        config_json = config_path+instance_type+str(num_r)+"_config.json"
        states = load_states(state_json)
        num_robots = len(states["robots"])
        traj = np.array([states["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
        collision_shape = load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2]
        FoV_beta = load_states(config_json)["fov_cbf_params"]["beta"] * np.pi/180
        FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]

        # Metric1: success rate
        success = instance_success(traj, collision_shape)
        print("success", success)
        # Metric2: avg number of neighbor in FoV
        avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta)
        mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
        print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)

# plots


