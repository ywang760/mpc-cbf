INPUT=(
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/circle.json"
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/2r/line.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/bend.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/circle.json"
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/cross_split.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line2.json"
     "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/line3.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/3r/triangle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/5r/circle.json"
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/5r/expand.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/6r/circle.json"
    #"/usr/src/mpc-cbf/workspace/experiments/config/baseline/6r/upward.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/8r/circle.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/baseline/8r/diverge.json"
)

DEFAULT_STATES_PATH="/usr/src/mpc-cbf/workspace/experiments/results/states.json"
BASE_CONFIG_FILE="/usr/src/mpc-cbf/workspace/experiments/config/base_config.json"
VIZ_OUTPUT_DIR="/usr/src/mpc-cbf/workspace/experiments/results/viz"
MAX_STEPS=5000

# Build the CBF examples once before running experiments
echo "Building CBF Formation Control example"
cd /usr/src/mpc-cbf/workspace/lib/cbf/build
make -j20 cbf_examples_CBFFormationControl_example

# Remove the existing states.json file if it exists
if [ -f ${DEFAULT_STATES_PATH} ]; then
    echo "Removing existing states.json file at ${DEFAULT_STATES_PATH}"
    rm ${DEFAULT_STATES_PATH}
fi

# Set logging level to debug
 export SPDLOG_LEVEL=info
 #export SPDLOG_LEVEL=debug
 #export SPDLOG_LEVEL=warn
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
    echo "Step 1: Running the CBF Formation Control example with configuration"
    ./cbf_examples_CBFFormationControl_example \
        --config_file ${config_file} \
        --write_filename ${DEFAULT_STATES_PATH} \
        --max_steps ${MAX_STEPS}

    # Step 2: Visualize the results
    echo "Step 2: Visualizing the results from the experiment"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/plot_results.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH} \
        --output_dir ${VIZ_OUTPUT_DIR}\
        --create_anim \
        --anim_format mp4
        # can add the --create_anim flag if needed

    # Step 3: Check for collisions and success
    echo "Step 3: Checking for collisions and success of the robot trajectories"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/metrics/collision_check.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH}
    
    echo "Completed processing: $(basename ${config_file})"
    echo ""
done

cd /usr/src/mpc-cbf/workspace/experiments/scripts
echo "All configurations processed successfully!"