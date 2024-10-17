import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

file = Path.home() / ("fovmpc/lib/particle_filter/build/samples.txt")
print("File path: ", file)
data = []
with open(str(file), 'r') as f:
    lines = f.readlines()

for l in lines:
    data.append(l)

N = 100
STEPS = len(lines) // N
covs = np.zeros((STEPS, N, 3))
for s in range(STEPS):
    for i in range(N):
        row = data[s*N + i].replace(',\n', '')
        row = tuple(map(float, row.split(',')))
        covs[s, i, :] = np.array(row)

fig, axs = plt.subplots(2, STEPS//2, figsize=(18,10))
row = 0
det = np.array([2.3, 1.5])
robot = np.array([1.0, 1.0])
for i in range(STEPS):
    if i >= STEPS//2:
        row = 1
    mean = np.mean(covs[i, :, :], axis=0)
    ax = axs[row, i-STEPS//2*row]
    ax.scatter(covs[i, :, 0], covs[i, :, 1], c='tab:blue')
    ax.scatter(mean[0], mean[1], c='tab:green', s=24, label="Estimate")
    ax.scatter(det[0], det[1], c='tab:red', s=24, label="Detection")
    ax.scatter(robot[0], robot[1], c='tab:orange', s=24, marker='X', label="Real")
    ax.set_xlim([-10.0, 15.0])
    ax.set_ylim([-10.0, 15.0])
    ax.legend()


plt.show()

fig, ax = plt.subplots(1, 1)
ax.scatter(covs[i, :, 0], covs[i, :, 1], c='tab:blue')
ax.scatter(mean[0], mean[1], c='tab:green', s=24, label="Estimate")
ax.scatter(det[0], det[1], c='tab:red', s=24, label="Detection")
ax.scatter(robot[0], robot[1], c='tab:orange', s=24, marker='X', label="Real")
ax.set_xlim([-10.0, 15.0])
ax.set_ylim([-10.0, 15.0])
ax.legend()
plt.show()


