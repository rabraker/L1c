FROM debian:bullseye

RUN apt-get update --yes \
    && apt-get install  --yes lsb-release software-properties-common wget gnupg2\
    && apt-get install  --yes lsb-release software-properties-common wget gnupg2 \
    && wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key|apt-key add - \
    && echo "deb http://apt.llvm.org/bullseye/ llvm-toolchain-bullseye-13 main" >> /etc/apt/sources.list \
    && apt-get update --yes \
    && apt-get install --yes \
    cmake \
    clang-format-13 \
    clang-tidy-13 \
    libopenblas-openmp-dev \
    python3 python3 python3-dev python3-pip \
    libfftw3-3 \
    libfftw3-dev \
    pkg-config \
    octave \
    octave-signal \
    liboctave-dev \
    git

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt && apt-get remove -y libopenblas0-pthread

COPY conanfile.txt conanfile.txt
COPY tools/conan_default_profile   /root/.conan/profiles/default
RUN  conan install .

# https://gitlab.com/rabraker/l1c/container_registry
# docker login registry.gitlab.com
# docker build -t registry.gitlab.com/rabraker/l1c .
# docker push registry.gitlab.com/rabraker/l1c
