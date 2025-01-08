cd ../../lib/mpc_cbf/cmake-build-release
cmake --build ./ --target mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example
date="01062025"
instances="circle"
max_r=10
max_exp=10
min_fov=120
max_fov=160
fov_step=20
for (( num_r = 2; num_r <= ${max_r}; num_r++ )); do
  for fov in $(seq ${min_fov} ${fov_step} ${max_fov}); do
    for (( exp_idx = 0; exp_idx < ${max_exp}; exp_idx++ )); do
        echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
        ./mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --fov ${fov} --write_filename "../../../experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_States_${exp_idx}.json"
    done
  done
done