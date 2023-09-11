FROM ubuntu:20.04
LABEL base_image="ubuntu:20.04"
LABEL version="2.0"
LABEL maintainer="Kevin Nathanael Ramanto"
LABEL maintainer.email="kevin@nalagenetics.com"

ENV TZ=Asia/Singapore
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone 

RUN export DEBIAN_FRONTEND=noninteractive \
    &&     apt-get update -y \
    &&     apt-get -y --no-install-recommends install git python-is-python3 python3-pip \
    &&     apt-get autoclean -y \
    &&     apt-get install -y apt-utils unzip wget build-essential zlib1g-dev libbz2-dev liblzma-dev libmath-libm-perl libpthread-stubs0-dev curl python3-pycurl git git-lfs gawk libncurses-dev libboost-all-dev libgtest-dev cmake libssl-dev libcurl3-dev\
    &&     pip install pandas numpy scipy argparse matplotlib 

RUN wget --max-redirect=1 https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2 \
    &&     tar -xf bcftools-1.15.1.tar.bz2 \
    &&     mkdir bcftools \
    &&     cd bcftools-1.15.1 \
    &&     ./configure --prefix=/bcftools \
    &&     make \
    &&     make install \
    &&     cd ../ \
    &&     rm bcftools-1.15.1.tar.bz2 \
    &&     apt-get install -y tabix \
    &&     apt-get install -y vim 


RUN mkdir GLIMPSE \
    &&     cd GLIMPSE \
    &&     wget --max-redirect=1 https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_chunk_static \
    &&     wget --max-redirect=1 https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_phase_static \
    &&     wget --max-redirect=1 https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_ligate_static \
    &&     wget --max-redirect=1 https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_split_reference_static \
    &&     wget --max-redirect=1 https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0/GLIMPSE2_concordance_static \
    &&     mv GLIMPSE2_chunk_static GLIMPSE2_chunk \
    &&     mv GLIMPSE2_phase_static GLIMPSE2_phase \
    &&     mv GLIMPSE2_ligate_static GLIMPSE2_ligate \
    &&     mv GLIMPSE2_split_reference_static GLIMPSE2_split_reference \
    &&     mv GLIMPSE2_concordance_static GLIMPSE2_concordance \
    &&     chmod +x GLIMPSE2_chunk \
    &&     chmod +x GLIMPSE2_phase \
    &&     chmod +x GLIMPSE2_ligate \
    &&     chmod +x GLIMPSE2_split_reference \
    &&     chmod +x GLIMPSE2_concordance \
    &&     cd ../

RUN mkdir wrapper

COPY ./GLIMPSE_phase.py /wrapper 

RUN pip install awscli 
ENV PATH="$PATH:/bcftools/bin"
ENV PATH="$PATH:/GLIMPSE/"
WORKDIR /data
