FROM ubuntu:22.04

ENV INSTALL_DIR="/opt"

RUN apt-get update; export DEBIAN_FRONTEND=noninteractive; apt-get -y --no-install-recommends install \
   wget \
   unzip \
   build-essential \
   gcc-multilib \
   netcat-openbsd && \
   export DEBIAN_FRONTEND=noninteractive; apt-get -y install openssh-client openssh-sftp-server openssh-server

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH
# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean --all -f -y

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/cnv_wes/bin:$PATH
USER 1000
