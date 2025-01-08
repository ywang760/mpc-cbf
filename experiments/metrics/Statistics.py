from Metrics import *
from ComputeCI import *
import colorsys

def generate_rgb_colors(num_colors):
    output = []
    num_colors += 1 # to avoid the first color
    for index in range(1, num_colors):
        incremented_value = 1.0 * index / num_colors
        output.append(colorsys.hsv_to_rgb(incremented_value, 0.5, 0.5))
    return np.asarray(output)

instance_type = "circle"
config_path = "../instances/"+instance_type+"_instances/"
result_path = "../instances/results/log"
date = "01062025"
min_r = 2
max_r = 10
min_fov=120
max_fov=160
fov_step=20
num_robot = range(min_r, max_r+1)
fovs = range(min_fov, max_fov+fov_step, fov_step)
num_exp = 10

success_rate_dict = {}
num_neighbor_in_fov_dict = {}
percent_neighbor_in_fov_dict = {}

success_rate_ci_dict={}
num_neighbor_in_fov_ci_dict={}
percent_neighbor_in_fov_ci_dict={}
for fov in fovs:
    success_rate_dict["fov"+str(fov)] = [[] for _ in num_robot]  # [entry,], averaged by M, M is the sample
    num_neighbor_in_fov_dict["fov"+str(fov)] = [[] for _ in num_robot]  # [entry, M]
    percent_neighbor_in_fov_dict["fov"+str(fov)] = [[] for _ in num_robot]  # [entry, M]

    success_rate_ci_dict["fov"+str(fov)] = {}
    num_neighbor_in_fov_ci_dict["fov"+str(fov)] = {}
    percent_neighbor_in_fov_ci_dict["fov"+str(fov)] = {}


for fov in fovs:
    for r_idx, num_r in enumerate(num_robot):
        print("r_idx: ", r_idx, "num_r: ", num_r)
        for exp_idx in range(num_exp):
            state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_States_"+str(exp_idx)+".json"
            config_json = config_path+instance_type+str(num_r)+"_fov"+str(fov)+"_config.json"
            states = load_states(state_json)
            num_robots = len(states["robots"])
            traj = np.array([states["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
            collision_shape = np.array(load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2])
            FoV_beta = load_states(config_json)["fov_cbf_params"]["beta"] * np.pi/180
            FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]

            # Metric1: success rate
            success = instance_success(traj, collision_shape)
            success_rate_dict["fov"+str(fov)][r_idx].append(float(success))

            # Metric2: avg number of neighbor in FoV
            num_neighbors = num_r - 1
            avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta)
            mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
            # print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)
            # print("mean_avg_num_neighbor_in_fov: ", mean_avg_num_neighbor_in_fov)
            num_neighbor_in_fov_dict["fov"+str(fov)][r_idx].append(mean_avg_num_neighbor_in_fov)
            percent_neighbor_in_fov_dict["fov"+str(fov)][r_idx].append(mean_avg_num_neighbor_in_fov / num_neighbors)

# plots
num_robot_arr = np.array(num_robot)
success_rate_mean_arr_list = []
success_rate_ci_arr_list = []
success_rate_label_list = []

num_neighbor_in_fov_mean_arr_list = []
num_neighbor_in_fov_ci_arr_list = []
num_neighbor_in_fov_label_list = []

percent_neighbor_in_fov_mean_arr_list = []
percent_neighbor_in_fov_ci_arr_list = []
percent_neighbor_in_fov_label_list = []

num_colors = 0
for fov in fovs:
    success_rate_ci_dict["fov"+str(fov)]["mean"], success_rate_ci_dict["fov"+str(fov)]["ci"] = CI_compute(np.array(success_rate_dict["fov"+str(fov)]))
    num_neighbor_in_fov_ci_dict["fov"+str(fov)]["mean"], num_neighbor_in_fov_ci_dict["fov"+str(fov)]["ci"] = CI_compute(np.array(num_neighbor_in_fov_dict["fov"+str(fov)]))
    percent_neighbor_in_fov_ci_dict["fov"+str(fov)]["mean"], percent_neighbor_in_fov_ci_dict["fov"+str(fov)]["ci"] = CI_compute(np.array(percent_neighbor_in_fov_dict["fov"+str(fov)]))

    success_rate_mean_arr_list.append(success_rate_ci_dict["fov"+str(fov)]["mean"])
    success_rate_ci_arr_list.append(success_rate_ci_dict["fov"+str(fov)]["ci"])
    success_rate_label_list.append("fov"+str(fov))
    num_neighbor_in_fov_mean_arr_list.append(num_neighbor_in_fov_ci_dict["fov"+str(fov)]["mean"])
    num_neighbor_in_fov_ci_arr_list.append(num_neighbor_in_fov_ci_dict["fov"+str(fov)]["ci"])
    num_neighbor_in_fov_label_list.append("fov"+str(fov))
    percent_neighbor_in_fov_mean_arr_list.append(percent_neighbor_in_fov_ci_dict["fov"+str(fov)]["mean"])
    percent_neighbor_in_fov_ci_arr_list.append(percent_neighbor_in_fov_ci_dict["fov"+str(fov)]["ci"])
    percent_neighbor_in_fov_label_list.append("fov"+str(fov))

    num_colors += 1

colors = generate_rgb_colors(num_colors)
assert len(colors) == len(success_rate_mean_arr_list)
list_CI_plot(num_robot_arr, success_rate_mean_arr_list, success_rate_ci_arr_list, success_rate_label_list, colors, xlabel="Number of Robots", ylabel="Success Rate", save_name="./multi_success_rate_ci.png")
list_CI_plot(num_robot_arr, num_neighbor_in_fov_mean_arr_list, num_neighbor_in_fov_ci_arr_list, num_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Number of Neighbors in FoV", save_name="./multi_nnif_ci.png")
list_CI_plot(num_robot_arr, percent_neighbor_in_fov_mean_arr_list, percent_neighbor_in_fov_ci_arr_list, percent_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Percentage of Neighbors in FoV", save_name="./multi_pnif_ci.png")
