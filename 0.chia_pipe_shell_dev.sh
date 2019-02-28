#!/bin/bash

# ChIA-PIPE
#         Starting point for launching ChIA-PIPE
# 2018
# The Jackson Laboratory for Genomic Medicine

## The help message:
function usage
{
    echo -e "usage: bash 0.chia_pipe_shell.sh -c CONF
    " 
}

## Parse the command-line argument (i.e., get the name of the config file)
while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


# Source the config file to get the parameter values
source ${conf}

## Add dependency dir to path
#workdir=$( pwd)
#dep_dir=$( cd ${dep_dir} && pwd )
#export PATH=${dep_dir}:${PATH}
#cd ${workdir}

# Set the output directory for writing files
out_dir="${run}"
mkdir -p ${out_dir}
cd ${out_dir}

# Print values
echo ${conf}
echo ${out_dir}

# Set the resource parameters for the computing cluster
# depending on the run type (miseq or hiseq)
n_thread=1
mem=8

## 1. Linker filtering
bash shell_scripts/1.filter_linker.sh --conf ${conf}

## 2a. Map "No linker" reads
bash shell_scripts/2.map.sh --conf ${conf} --tag_name none

## 2b. Map "Linker, one tag" reads
bash shell_scripts/2.map.sh --conf ${conf} --tag_name singlelinker.single

## 2c. Map "Linker, two tag" reads
bash shell_scripts/2.map.sh --conf ${conf} --tag_name singlelinker.paired

## 3. Peak calling
bash shell_scripts/3.call_peaks.sh --conf ${conf}

## 4. QC table generation
bash shell_scripts/4.extract_summary_stats.sh --conf ${conf}

## 5. Allele-specific analysis
bash shell_scripts/5.phase_loops.sh --conf ${conf}

