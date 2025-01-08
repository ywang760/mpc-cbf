cd ../../tools/
date="01062025"
instances="circle"
min_r=2
max_r=10
min_fov=120
max_fov=160
fov_step=20
exp_id=0
for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
  for fov in $(seq ${min_fov} ${fov_step} ${max_fov}); do
    echo Animating simulation of ${num_r} robots...
    python3 Plotter.py --config_filename "../experiments/instances/${instances}_instances/${instances}${num_r}_fov${fov}_config.json" --states_filename "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_States_${exp_id}.json" --output_filename "../experiments/instances/results/log${date}/${instances}${num_r}.mp4"
  done
done