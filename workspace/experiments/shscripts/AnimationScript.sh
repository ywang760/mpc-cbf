#cd ../../tools/
#date="01142025"
#instances="circle"
##instances="formation"
#min_r=2
#max_r=10
#min_fov=120
#max_fov=120
#fov_step=120
#min_slack_decay=0.2
#max_slack_decay=0.3
#slack_decay_step=0.1
#default_slack_decay=0.1
#exp_id=0
#for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
#  for fov in $(seq ${min_fov} ${fov_step} ${max_fov}); do
#    echo Animating simulation of ${num_r} robots...
#    python3 Plotter.py --fov ${fov} --config_filename "../experiments/instances/${instances}_instances/${instances}${num_r}_config.json" --states_filename "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${default_slack_decay}_States_${exp_id}.json" --output_video "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${default_slack_decay}.mp4" --output_figure "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${default_slack_decay}.png"
#  done
#done

#for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
#  for fov in 120 140; do
#    for slack_decay in $(seq ${min_slack_decay} ${slack_decay_step} ${max_slack_decay}); do
#      echo Animating simulation of ${num_r} robots...
#      python3 Plotter.py --config_filename "../experiments/instances/${instances}_instances/${instances}${num_r}_config.json" --states_filename "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${slack_decay}_States_${exp_id}.json" --output_filename "../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${slack_decay}.mp4"
#    done
#  done
#done