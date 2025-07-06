INPUT=(
    # "/usr/src/mpc-cbf/workspace/experiments/config/formation/Separate2.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/formation/Cross3.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/formation/DiagonalSwitch3.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/formation/Horizontal3.json"
    # "/usr/src/mpc-cbf/workspace/experiments/config/formation/Circle5.json"
    "/usr/src/mpc-cbf/workspace/experiments/config/formation/Circle6.json"
)

DEFAULT_STATES_PATH="/usr/src/mpc-cbf/workspace/experiments/results/formation/states.json"
MAX_STEPS=1500
BASE_CONFIG_FILE="/usr/src/mpc-cbf/workspace/experiments/config/formation/base_config.json"

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

# Iterate over each configuration file
for config_file in "${INPUT[@]}"; do
    echo "=========================================="
    echo "Processing configuration: $(basename ${config_file})"
    echo "=========================================="
    
    # Step 0: Visualize the config
    echo "Step 0: Processing and visualizing the configuration file"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/preprocess.py \
        --base_config_file ${BASE_CONFIG_FILE} \
        --task_config_file ${config_file}

    # Step 1: Run the experiment
    echo "Step 1: Running the CBF Formation Control example with configuration"
    ./cbf_examples_CBFFormationControl_example \
        --config_file ${config_file} \
        --max_steps ${MAX_STEPS}

    # Step 2: Visualize the results
    echo "Step 2: Visualizing the results from the experiment"
    python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/plot_results.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH} \
        --output_dir /usr/src/mpc-cbf/workspace/experiments/results/formation/viz

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