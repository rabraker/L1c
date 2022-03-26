FROM gcc

RUN apt-get update --yes \
    && apt-get install  --yes lsb-release software-properties-common \
    && apt-get update --yes \
    && apt-get install --yes \
    cmake         \
    pkg-config \
    libopenblas-dev \
    python3 python3-numpy python3-scipy python3-dev python3-setuptools \
    libsubunit0 \
    libfftw3-3 libfftw3-dev \
    pkg-config

RUN apt-get install -y check

RUN git clone https://github.com/DaveGamble/cJSON.git && cd cJSON \
    && mkdir build && cd build && cmake .. && make -j$(nproc) \
    && make install \
    && ls /usr/local/include \
    && ls /usr/local/lib \
    && ldconfig

# https://gitlab.com/rabraker/l1c/container_registry
# docker login registry.gitlab.com
# docker build -t registry.gitlab.com/rabraker/l1c .
# docker push registry.gitlab.com/rabraker/l1c
