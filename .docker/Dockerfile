FROM ubuntu:22.04

# Set the working directory inside the container
WORKDIR /usr/src/mpc-cbf

# Set environment variables to make installation non-interactive
ENV DEBIAN_FRONTEND=noninteractive
# Install build tools
RUN apt-get update && apt-get install -y \
    python3 \
    python-is-python3 \
    python3-pip \
    python3.10-venv \
    libboost-all-dev \
    libeigen3-dev \
    libginac-dev \
    protobuf-compiler \
    lsb-release \
    curl \
    build-essential \
    cmake \
    doxygen \
    git \
    ninja-build \
    ffmpeg

# Copy the current directory contents into the container
COPY . .

# Initialize the submodules
RUN git submodule update --init --recursive

# Configure Git to trust the workspace path
RUN git config --global --add safe.directory /usr/src/mpc-cbf

# Setup CPLEX
ENV CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio2211
ENV PATH=$CPLEX_HOME/cplex/bin/x86-64_linux:$PATH
ENV CPLEX_INSTALLER=cplex_studio*.linux*.bin
RUN echo "Starting CPLEX installation"
WORKDIR /usr/src/mpc-cbf/.docker
RUN set -x && \
    if [ -f cplex/$CPLEX_INSTALLER ]; then \
    echo "Installing CPLEX" && \
    chmod a+rwx cplex/$CPLEX_INSTALLER && \
    cplex/$CPLEX_INSTALLER -f cplex.installer.properties; \
    else \
    echo "No installer found for CPLEX"; \
    fi
RUN rm -rf /usr/src/mpc-cbf/.docker/cplex

# Initialize shell
WORKDIR /usr/src/mpc-cbf
ENV DEBIAN_FRONTEND=dialog
CMD ["bash"]