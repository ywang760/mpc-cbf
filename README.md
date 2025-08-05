# MPC-CBF: Model Predictive Control with Control Barrier Functions

## Getting Started

Clone the repository:

```bash
git clone https://github.com/ywang760/mpc-cbf.git
```

### Installation

The most convenient way to use the code is through a [docker container](https://www.docker.com/).

```bash
docker pull wangyt163/mpc-cbf:latest
```

The above container includes all the necessary dependencies to run the code. If you would like to build the container yourself, see [here](.docker/README.md).

If you'd like to run the code natively, key dependencies include:
- Eigen3
- Boost
- Ginac
- ILOG CPLEX Optimization Studio

### Running the code

In the root directory, run the following command:

```bash
docker compose -f .docker/docker-compose.yml up
```

This would start a container named `dev-container` and keep it running in the background.

Then, click on the icon in the bottom left corner of VSCode and select `Attach to Running Container...`, select the `mpc-cbf-dev-container` container.
You should be able to see the code in the `/usr/src/mpc-cbf` directory and run the code from there.
The `workspace` directory is mounted in the container, so any changes you make in the directory will be persisted throughout sessions.

(Make sure to have the [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed in VSCode.)

In the container, you can use the following commands to compile a library (take the CBF library as an example):

```bash
cd /usr/src/mpc-cbf/workspace/lib/cbf
mkdir build
cd build
cmake .. -G Ninja (or use `cmake ..` to use the default generator)
ninja (or `make` to use the default generator)
```

This would compile the CBF library with ninja (installed in the container by default).

## Code Structure

The code is organized as follows:


## Experiments

### Baseline experiments using CBF

- 