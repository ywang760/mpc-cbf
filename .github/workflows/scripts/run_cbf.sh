INPUT=(
    "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/2r/circle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/2r/line.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/bend.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/circle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/cross_split.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/line.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/line2.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/line3.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/3r/triangle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/5r/circle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/5r/expand.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/6r/circle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/6r/upward.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/8r/circle.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/8r/circle2.json"
    # "$GITHUB_WORKSPACE/workspace/experiments/config/baseline/8r/diverge.json"
)

DEFAULT_STATES_PATH="$GITHUB_WORKSPACE/workspace/experiments/results/states.json"
BASE_CONFIG_FILE="$GITHUB_WORKSPACE/workspace/experiments/config/base_config.json"
VIZ_OUTPUT_DIR="$GITHUB_WORKSPACE/workspace/experiments/results/ci_plots"
MAX_STEPS=300

# Build the CBF examples once before running experiments
echo "Building CBF Formation Control example"
cd $GITHUB_WORKSPACE/workspace/lib/cbf/build
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
    python3 $GITHUB_WORKSPACE/workspace/experiments/python/preprocess.py \
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
    python3 $GITHUB_WORKSPACE/workspace/experiments/python/visualization/plot_results.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH} \
        --output_dir ${VIZ_OUTPUT_DIR}\
        --create_anim \
        --anim_format mp4
        # can add the --create_anim flag if needed

    # Step 3: Check for collisions and success
    echo "Step 3: Checking for collisions and success of the robot trajectories"
    python3 $GITHUB_WORKSPACE/workspace/experiments/python/metrics/collision_check.py \
        --config ${config_file} \
        --states ${DEFAULT_STATES_PATH}
    
    echo "Completed processing: $(basename ${config_file})"
    echo ""
done

echo "All configurations processed successfully!"