from Metrics import *
from ComputeCI import *

result_path = "../instances/results/log"
config_path = "../instances/"
date = "01042025"
instance_type = "circle"
min_r = 2
max_r = 10
num_robot = range(min_r, max_r+1)
num_exp = 5

success_rate = [[] for _ in num_robot]  # [entry,], averaged by M, M is the sample
num_neighbor_in_fov = [[] for _ in num_robot]  # [entry, M]
percent_neighbor_in_fov = [[] for _ in num_robot]  # [entry, M]

for r_idx, num_r in enumerate(num_robot):
    print("r_idx: ", r_idx, "num_r: ", num_r)
    for exp_idx in range(num_exp):
        state_json = result_path+date+"/"+instance_type+str(num_r)+"States_"+str(exp_idx)+".json"
        config_json = config_path+instance_type+str(num_r)+"_config.json"
        states = load_states(state_json)
        num_robots = len(states["robots"])
        traj = np.array([states["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
        collision_shape = np.array(load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2])
        FoV_beta = load_states(config_json)["fov_cbf_params"]["beta"] * np.pi/180
        FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]

        # Metric1: success rate
        success = instance_success(traj, collision_shape)
        success_rate[r_idx].append(float(success))

        # Metric2: avg number of neighbor in FoV
        num_neighbors = num_r - 1
        avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta)
        mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
        # print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)
        # print("mean_avg_num_neighbor_in_fov: ", mean_avg_num_neighbor_in_fov)
        num_neighbor_in_fov[r_idx].append(mean_avg_num_neighbor_in_fov)
        percent_neighbor_in_fov[r_idx].append(mean_avg_num_neighbor_in_fov / num_neighbors)

# plots
print("success_rate: ", success_rate)
print("num_neighbor_in_fov: ", num_neighbor_in_fov)
print("percent_neighbor_in_fov: ", percent_neighbor_in_fov)

success_rate_arr = np.array(success_rate)
num_neighbor_in_fov_arr = np.array(num_neighbor_in_fov)
percent_neighbor_in_fov_arr = np.array(percent_neighbor_in_fov)

num_robot_arr = np.array(num_robot)
success_mean, success_ci = CI_compute(success_rate_arr)
nnif_mean, nnif_ci = CI_compute(num_neighbor_in_fov_arr)
pnif_mean, pnif_ci = CI_compute(percent_neighbor_in_fov_arr)

CI_plot(num_robot_arr, success_mean, success_ci, save_name="./success_rate_ci.png", xlabel="#robots", ylabel="Success Rate")
CI_plot(num_robot_arr, nnif_mean, nnif_ci, save_name="./nnif_ci.png", xlabel="#robots", ylabel="#neighbors in FoV")
CI_plot(num_robot_arr, pnif_mean, pnif_ci, save_name="./pnif_ci.png", xlabel="#robots", ylabel="%neighbors in FoV")