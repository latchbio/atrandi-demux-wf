from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

run curl -L -O \
    https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Mambaforge-Linux-x86_64.sh -b \
    && rm -f Mambaforge-Linux-x86_64.sh

run /root/mambaforge/bin/mamba create -n snakemake python=3.11
env PATH /root/mambaforge/envs/snakemake/bin:$PATH

run pip install snakemake snakemake_executor_plugin_latch==0.1.2 snakemake_storage_plugin_latch==0.1.9

# Latch SDK
# DO NOT REMOVE
run pip install latch==2.54.0.a9
run mkdir /opt/latch

copy wf /root/wf
copy latch_metadata /root/latch_metadata
copy Dockerfile /root/Dockerfile
copy Snakefile.demux /root/Snakefile.demux
copy version /root/version
# copy scripts /root/scripts
copy config.yaml /root/config.yaml
copy workflow /root/workflow

arg tag
env FLYTE_INTERNAL_IMAGE $tag
workdir /root
