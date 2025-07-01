INPUT=/usr/src/mpc-cbf/workspace/experiments/config/formation/robots2_2.json

# Step 0: Visualize the config
echo "Visualizing the configuration file"
python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/plot_config.py \
    ${INPUT}

# Step 1: Run the experiment
echo "Running the CBF Formation Control example with configuration"
cd /usr/src/mpc-cbf/workspace/lib/cbf/build
make -j4 cbf_examples_CBFFormationControl_example
./cbf_examples_CBFFormationControl_example \
    --config_file ${INPUT}

# Step 2: Visualize the results
echo "Visualizing the results from the experiment"
python3 /usr/src/mpc-cbf/workspace/experiments/python/visualization/plot_results.py \
    --config ${INPUT} \
    --states /usr/src/mpc-cbf/workspace/experiments/results/formation/states.json \
    --output_dir /usr/src/mpc-cbf/workspace/experiments/results/formation/viz