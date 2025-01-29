from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

run curl -L -O \
    https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Mambaforge-Linux-x86_64.sh -b \
    && rm -f Mambaforge-Linux-x86_64.sh

run /root/mambaforge/bin/mamba create -n snakemake python=3.11

env PATH /root/mambaforge/bin:$PATH
env PATH /root/mambaforge/envs/snakemake/bin:$PATH

copy envs /root/envs

run mamba env create -n minimap -f envs/minimap.yaml
run mamba env create -n pandas -f envs/pandas.yaml
run mamba env create -n pheniqs -f envs/pheniqs.yaml
run mamba env create -n plot-observed -f envs/plot-observed.yaml
run mamba env create -n seqtk -f envs/seqtk.yaml

# Latch SDK
# DO NOT REMOVE
run pip install "latch[snakemake]"==2.55.0a6
run mkdir /opt/latch

copy wf /root/wf
copy latch_metadata /root/latch_metadata
copy Dockerfile /root/Dockerfile
copy Snakefile.demux /root/Snakefile.demux
copy version /root/version
copy config.yaml /root/config.yaml
copy workflow /root/workflow

arg tag
env FLYTE_INTERNAL_IMAGE $tag
workdir /root
