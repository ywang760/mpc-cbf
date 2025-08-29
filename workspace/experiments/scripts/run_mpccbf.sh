INPUT=(
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/line.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/bend.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/cross_split.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line2.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line3.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/triangle.json"
    "/usr/src/mpc-cbf/workspace/experiments/config/baseline/5r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/5r/expand.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/6r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/6r/upward.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/8r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/8r/diverge.json"
)

DEFAULT_STATES_PATH="/usr/src/mpc-cbf/workspace/experiments/results/states.json"
BASE_CONFIG_FILE="/usr/src/mpc-cbf/workspace/experiments/config/base_config.json"
VIZ_OUTPUT_DIR="/usr/src/mpc-cbf/workspace/experiments/results/mpccbf_viz"
SIM_RUNTIME=10.0

# Build the MPC CBF examples once before running experiments
echo "Building MPC CBF Formation Control example (Make)"
SRC_DIR="/usr/src/mpc-cbf/workspace/lib/mpc_cbf"
BUILD_DIR="${SRC_DIR}/build"
# 如果 build 目录不存在就先 cmake
if [ ! -d "$BUILD_DIR" ]; then
  cmake -S "$SRC_DIR" -B "$BUILD_DIR" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
fi
# 增量编译
make -C "$BUILD_DIR" -j"$(nproc)"

# Remove the existing states.json file if it exists
if [ -f ${DEFAULT_STATES_PATH} ]; then
    echo "Removing existing states.json file at ${DEFAULT_STATES_PATH}"
    rm ${DEFAULT_STATES_PATH}
fi

# Set logging level to debug
export SPDLOG_LEVEL=info

# Iterate over each configuration file
for config_file in "${INPUT[@]}"; do
    echo "=========================================="
    echo "Processing configuration: $(basename ${config_file})"
    echo "=========================================="
    
    # Step 0: Preprocesss
    echo "Step 0: Processing the configuration file"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/preprocess.py \
        --base_config_file ${BASE_CONFIG_FILE} \
        --task_config_file ${config_file}

    # Step 1: Run the experiment
    echo "Step 1: Running the MPC CBF Formation Control example with configuration"
    /usr/src/mpc-cbf/workspace/lib/mpc_cbf/build/mpc_cbf_examples_MPCCBFFormationControl_example \
        --config_file "${config_file}" \
        --write_filename "${DEFAULT_STATES_PATH}" \
        --sim_runtime "${SIM_RUNTIME}"

    # Step 3: Check for collisions and success
    echo "Step 3: Checking for collisions and success of the robot trajectories"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/metrics/collision_check.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH}
    
    echo "Completed processing: $(basename ${config_file})"
    echo ""
    
    # Step 2: Visualize the results
    echo "Step 2: Visualizing the results from the experiment"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/plot_results.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH} \
        --output_dir ${VIZ_OUTPUT_DIR} \
        --create_anim
        # can add the --create_anim flag if needed
done

cd /usr/src/mpc-cbf/workspace/experiments/scripts
echo "All configurations processed successfully!"