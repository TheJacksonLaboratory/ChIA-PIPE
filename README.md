# README
# ChIA-PIPE
## The Jackson Laboratory for Genomic Medicine
## 26 November 2018


## Download ChIA-PIPE
git clone git@github.com:TheJacksonLaboratory/chia_pipe.git

## Or download directly:
wget https://github.com/TheJacksonLaboratory/chia_pipe.zip

## Install ChIA-PIPE dependencies
dep_dir="/path/to/local/install/dependencies"


bash local_install_chia_pipe_dependencies.sh -i ${dep_dir}

## Download test data from Dropbox
https://www.dropbox.com/sh/kiunzdfj74dpamh/AAB3u4vVGHIUSZoiPSjEYuYva?dl=0

mkdir -p LDK0004-ds

cp LDK0004-ds_*.fastq.gz LDK0004-ds

## Launch ChIA-PIPE
qsub -F "--conf my_config_file.sh" 0.chia_pipe.pbs
