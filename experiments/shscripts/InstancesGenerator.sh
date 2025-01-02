cd ../instances/
min_r=2
max_r=8
for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
  python3 CircleInstances.py --num_robots ${num_r}
done