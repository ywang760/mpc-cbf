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
    ninja-build
    
# Setup g++-9 as the compiler
# RUN apt-get update && apt-get install -y \
#     g++-9 \
#     gcc-9
# RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100 && \
#     update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100
# ENV CC=gcc-9
# ENV CXX=g++-9

# Copy the current directory contents into the container
COPY . .

# Initialize the submodules
RUN git submodule update --init --recursive

# Setup CPLEX
ENV CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio2211
ENV PATH=$CPLEX_HOME/cplex/bin/x86-64_linux:$PATH
ENV LD_LIBRARY_PATH=$CPLEX_HOME/cplex/lib/x86-64_linux:${LD_LIBRARY_PATH:-""}
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

# # Build the cbf library
# WORKDIR /usr/src/mpc-cbf/lib/cbf
# RUN mkdir build && cd build && cmake -G Ninja .. && ninja

# Initialize shell
WORKDIR /usr/src/mpc-cbf
ENV DEBIAN_FRONTEND=dialog
CMD ["bash"]