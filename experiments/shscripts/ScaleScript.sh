cd ../../lib/mpc_cbf/cmake-build-release
cmake --build ./ --target mpc_cbf_examples_BezierIMPCCBFXYYaw_example
date="01012025"
instances="circle"
max_r=8
max_exp=3
for (( num_r = 2; num_r <= ${max_r}; num_r++ )); do
  for (( exp_idx = 0; exp_idx < ${max_exp}; exp_idx++ )); do
      echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
      ./mpc_cbf_examples_BezierIMPCCBFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --write_filename "../../../experiments/instances/results/log${date}/${instances}${num_r}States_${exp_idx}.json"
  done
done