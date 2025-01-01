cd ../../tools/
date="01012025"
instances="circle"
max_r=8
exp_id=0
for (( num_r = 2; num_r <= ${max_r}; num_r++ )); do
  echo Animating simulation of ${num_r} robots...
  python3 Plotter.py --config_filename "../experiments/instances/${instances}${num_r}_config.json" --states_filename "../experiments/instances/results/log${date}/${instances}${num_r}States_${exp_id}.json" --output_filename "../experiments/instances/results/log${date}/${instances}${num_r}.mp4"
done