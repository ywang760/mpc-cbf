cd ../instances/
min_r=2
max_r=10
min_fov=120
max_fov=180
fov_step=20
for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
  for fov in $(seq ${min_fov} ${fov_step} ${max_fov}); do
    python3 CircleInstances.py --num_robots ${num_r} --fov ${fov}
  done
done