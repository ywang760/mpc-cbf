from Metrics import *
from ComputeCI import *
import colorsys

def generate_rgb_colors(num_colors):
    output = []
    num_colors += 1 # to avoid the first color
    for index in range(1, num_colors):
        incremented_value = 1.0 * index / num_colors
        output.append(colorsys.hsv_to_rgb(incremented_value, 0.5, 0.7))
    return np.asarray(output)

instance_type = "circle"
# instance_type = "formation"
config_path = "../instances/"+instance_type+"_instances/"
result_path = "/media/lishuo/ssd/RSS2025_results/log"
date = "01102025"
num_exp = 10
min_r = 2
max_r = 10
min_fov=120
max_fov=360
fov_step=60
min_slack_decay=0.2
max_slack_decay=0.3
slack_decay_step=0.1
default_slack_decay=0.2
decay_exp_fovs = [120, 180, 240]
num_robot = range(min_r, max_r+1)
fovs = range(min_fov, max_fov+fov_step, fov_step)
slack_decays = [0.1, 0.3, 0.4]
print("slack decays: ", slack_decays)

experiment_key=[]
experiment_slack_decay_key=[]
experiment_fov_key=[]
for fov in fovs:
    experiment_key.append("fov"+str(fov)+"decay"+str(default_slack_decay))
    experiment_slack_decay_key.append(default_slack_decay)
    experiment_fov_key.append(fov)
for fov in decay_exp_fovs:
    for decay in slack_decays:
        experiment_key.append("fov"+str(fov)+"decay"+str(decay))
        experiment_slack_decay_key.append(decay)
        experiment_fov_key.append(fov)

success_rate_dict = {}
num_neighbor_in_fov_dict = {}
percent_neighbor_in_fov_dict = {}

success_rate_ci_dict={}
num_neighbor_in_fov_ci_dict={}
percent_neighbor_in_fov_ci_dict={}
for key in experiment_key:
    success_rate_dict[key] = [[] for _ in num_robot]  # [entry,], averaged by M, M is the sample
    num_neighbor_in_fov_dict[key] = [[] for _ in num_robot]  # [entry, M]
    percent_neighbor_in_fov_dict[key] = [[] for _ in num_robot]  # [entry, M]

    success_rate_ci_dict[key] = {}
    num_neighbor_in_fov_ci_dict[key] = {}
    percent_neighbor_in_fov_ci_dict[key] = {}

for i in range(len(experiment_key)):
    key = experiment_key[i]
    fov = experiment_fov_key[i]
    slack_decay = experiment_slack_decay_key[i]
    for r_idx, num_r in enumerate(num_robot):
        print("r_idx: ", r_idx, "num_r: ", num_r)
        for exp_idx in range(num_exp):
            state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_decay"+str(slack_decay)+"_States_"+str(exp_idx)+".json"
            config_json = config_path+instance_type+str(num_r)+"_config.json"
            states = load_states(state_json)
            goals = np.array(load_states(config_json)["tasks"]["sf"])
            num_robots = len(states["robots"])
            traj = np.array([states["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
            collision_shape = np.array(load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2])
            FoV_beta = fov * np.pi/180
            FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]

            # Metric1: success rate
            goal_radius = 1.6
            success = instance_success(traj, goals, goal_radius, collision_shape)
            print("success: ", success)
            success_rate_dict[key][r_idx].append(float(success))

            # Metric2: avg number of neighbor in FoV
            num_neighbors = num_r - 1
            avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta)
            mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
            # print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)
            # print("mean_avg_num_neighbor_in_fov: ", mean_avg_num_neighbor_in_fov)
            num_neighbor_in_fov_dict[key][r_idx].append(mean_avg_num_neighbor_in_fov)
            percent_neighbor_in_fov_dict[key][r_idx].append(mean_avg_num_neighbor_in_fov / num_neighbors)

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
for key in experiment_key:
    success_rate_ci_dict[key]["mean"], success_rate_ci_dict[key]["ci"] = CI_compute(np.array(success_rate_dict[key]))
    num_neighbor_in_fov_ci_dict[key]["mean"], num_neighbor_in_fov_ci_dict[key]["ci"] = CI_compute(np.array(num_neighbor_in_fov_dict[key]))
    percent_neighbor_in_fov_ci_dict[key]["mean"], percent_neighbor_in_fov_ci_dict[key]["ci"] = CI_compute(np.array(percent_neighbor_in_fov_dict[key]))

    success_rate_mean_arr_list.append(success_rate_ci_dict[key]["mean"])
    success_rate_ci_arr_list.append(success_rate_ci_dict[key]["ci"])
    success_rate_label_list.append(key)
    num_neighbor_in_fov_mean_arr_list.append(num_neighbor_in_fov_ci_dict[key]["mean"])
    num_neighbor_in_fov_ci_arr_list.append(num_neighbor_in_fov_ci_dict[key]["ci"])
    num_neighbor_in_fov_label_list.append(key)
    percent_neighbor_in_fov_mean_arr_list.append(percent_neighbor_in_fov_ci_dict[key]["mean"])
    percent_neighbor_in_fov_ci_arr_list.append(percent_neighbor_in_fov_ci_dict[key]["ci"])
    percent_neighbor_in_fov_label_list.append(key)

    num_colors += 1

colors = generate_rgb_colors(num_colors)
assert len(colors) == len(success_rate_mean_arr_list)
list_CI_plot(num_robot_arr, success_rate_mean_arr_list, success_rate_ci_arr_list, success_rate_label_list, colors, xlabel="Number of Robots", ylabel="Success Rate", save_name="./multi_success_rate_ci.png")
list_CI_plot(num_robot_arr, num_neighbor_in_fov_mean_arr_list, num_neighbor_in_fov_ci_arr_list, num_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Number of Neighbors in FoV", save_name="./multi_nnif_ci.png")
list_CI_plot(num_robot_arr, percent_neighbor_in_fov_mean_arr_list, percent_neighbor_in_fov_ci_arr_list, percent_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Percentage of Neighbors in FoV", save_name="./multi_pnif_ci.png")
