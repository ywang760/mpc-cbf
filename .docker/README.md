# Docker Usage

## Obtaining the container

TODO: include instructions on `docker pull`

## Building the container

A `.docker` directory is included in the root of the repository. This directory contains the `Dockerfile` and any other necessary files for building the container. The `Dockerfile` specifies the base image and the dependencies required to run the code.

To build the container, run the following command from the root of the repository:

```bash
docker build -f .docker/Dockerfile -t mpc-cbf .
```

This command will build the container and tag it with the name `mpc-cbf`.

This explicit build process can be skipped if using docker-compose. See [Running the container using docker-compose](#running-the-container-using-docker-compose) for more information.

### Including ILOG CPLEX Optimization Studio

To use the ILOG CPLEX Optimization Studio, you must have a valid license file. See the [IBM website](https://www.ibm.com/products/ilog-cplex-optimization-studio) for more information.

Once a license is obtained:

1. Download the installer `cplex_studioXXXX.linux-x86-64.bin` from the IBM website, replacing `XXXX` with the version number (e.g. 2211).
2. Place the installer in the `.docker/cplex` directory.

The directory should look like this:

```
.docker
├── cplex
│   └── cplex_studioXXXX.linux-x86-64.bin
|   └── cplex.installer.properties
```

## Using the container

The container can be run in either interactive mode (through `docker run -it`) or using docker-compose. The following sections provide instructions for both methods.

### Running the container in interactive mode (`docker run`)

To run the container, use the following command:

```bash
docker run -it mpc-cbf
```

To break down the command:

- `-it` runs the container in interactive mode. `mpc-cbf` is the name of the image.

Optional flags:

- `--name=<container-name>` assigns a name to the container.

This command will start the container and open a shell prompt. You can now run the code from the container.
However, any data created in the container will be lost when the container is stopped. To persist data, mount additional volumes.

**Note that this command will not build the binaries in the container. To run the scripts, you must manually build using the standard `cmake .. && make` procedure in the `build` directory.**

### Running the container using docker-compose

There is a `docker-compose.yml` file in the `.docker` directory that can be used to run the container. It contains a builder service and multiple runner services. The builder service builds the binaries in `lib/examples`, and shares the build directory with the runner services using a `build` volume. The runner services require the builder service to be launched first.

For more information on docker-compose, see the [official documentation](https://docs.docker.com/compose/).

#### To launch the builder service:

```bash
docker compose -f .docker/docker-compose.yml up builder
```

An optional `--build` flag can be added to build the container (skipping [explicity build step](#building-the-container)) or to rebuild the container.

```bash
docker compose -f .docker/docker-compose.yml up --build builder
```

#### To launch a runner service:

```bash
docker compose -f .docker/docker-compose.yml up <service-name>
```
