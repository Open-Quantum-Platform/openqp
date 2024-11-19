# Use a base Linux image with Python support
FROM ubuntu:20.04

# Set timezone to Asia/Seoul to avoid timezone prompt
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Asia/Seoul /etc/localtime && \
    echo "Asia/Seoul" > /etc/timezone

# Install dependencies, including LAPACK and BLAS, with --fix-missing
RUN apt-get update && apt-get install -y --fix-missing \
    gcc g++ gfortran cmake ninja-build \
    openmpi-bin libopenmpi-dev \
    python3-pip wget git \
    libblas-dev liblapack-dev \
    && apt-get clean

# Install Python dependencies
RUN pip3 install cffi

# Install CMake 3.25.2
WORKDIR /tmp
RUN wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz
RUN tar -zxvf cmake-3.25.2-linux-x86_64.tar.gz
RUN mv cmake-3.25.2-linux-x86_64 /opt/cmake
ENV PATH="/opt/cmake/bin:$PATH"

# Download and compile OpenQP
RUN git clone https://github.com/Open-Quantum-Platform/openqp.git /opt/openqp
WORKDIR /opt/openqp
RUN cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF
RUN ninja -C build install
RUN cd pyoqp && pip3 install .

# Install dftd4
WORKDIR /opt
RUN wget https://github.com/dftd4/dftd4/releases/download/v3.3.0/dftd4-3.3.0-linux-x86_64.tar.xz
RUN tar -xf dftd4-3.3.0-linux-x86_64.tar.xz && rm dftd4-3.3.0-linux-x86_64.tar.xz
RUN mv dftd4-3.3.0 /opt/dftd4
ENV PATH="/opt/dftd4/bin:$PATH"
ENV DFTD4_BIN=/opt/dftd4/bin/dftd4

# Install dftd4 via pip
RUN pip3 install dftd4

# Set environment variables, initializing LD_LIBRARY_PATH if not defined
ENV OPENQP_ROOT=/opt/openqp
ENV OMP_NUM_THREADS=4

# Run tests to confirm installation
RUN openqp --run_tests all

# Set entrypoint if required
ENTRYPOINT ["bash"]
