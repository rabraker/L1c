FROM gcc

RUN apt-get update --yes \
    && apt-get install  --yes lsb-release software-properties-common \
    && apt-get update --yes \
    && apt-get install --yes \
    cmake         \
    libopenblas-openmp-dev \
    python3 python3 python3-dev python3-pip \
    libfftw3-3 \
    libfftw3-dev \
    pkg-config

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY conanfile.txt conanfile.txt
COPY tools/conan_default_profile   /root/.conan/profiles/default
RUN  conan install .

# https://gitlab.com/rabraker/l1c/container_registry
# docker login registry.gitlab.com
# docker build -t registry.gitlab.com/rabraker/l1c .
# docker push registry.gitlab.com/rabraker/l1c
