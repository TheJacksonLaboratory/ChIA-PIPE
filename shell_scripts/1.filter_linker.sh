#!/bin/bash

# ChIA-PIPE
#         Step 1: Filter the linker sequence
#         (Requires single-linker ChIA-PET data)
# 2018
# The Jackson Laboratory for Genomic Medicine


## The help message:
function usage
{
    echo -e "usage: bash 1.filter_linker.sh --conf ${conf} 
    "
}

## Parse arguments from the command line
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


## Source config file
source ${conf}

# Add dependency dir to path
export PATH=${dep_dir}:${PATH}

## Create name of the log file
log_file=1.${run}.filter_linker.log

## Report arguments parsed
echo "
Arguments:
    main_prog=${main_prog}
    run=${run}
    data_dir=${data_dir}
    out_dir=${out_dir}
" >> ${log_file}


## Perform linker detection and generate different categories of fastq files
# Report linker detection start
echo "
`date` --- Linker detection started on: ---
    ${data_dir}/${r1_fastq},
    ${data_dir}/${r2_fastq}
" >> ${log_file}


## Linker filtering
if [ ${experiment_type} == 'ChIA-PET' ]; then
    # Linker filtering for ChIA-PET experiment

    if [ ${linker_a} == "none" ]; then
        # Default linker sequence
        # Single linker with default sequence
        ${main_prog} stag -W -T ${min_tag_len} -t ${n_thread} -O ${run} \
            ${data_dir}/${r1_fastq} ${data_dir}/${r2_fastq} \
            2>> ${log_file}
    
    elif [ ${linker_a} != "none" ]; then
        # Custom linker sequence
        if [ ${linker_b} == "none" ]; then
            # Single linker with custom sequence
            ${main_prog} stag -W -T ${min_tag_len} -t ${n_thread} -O ${run} \
                -A ${linker_a} \
                ${data_dir}/${r1_fastq} ${data_dir}/${r2_fastq} \
                2>> ${log_file}
        elif [ ${linker_b} != "none" ]; then
            # Two linkers with custom sequence
                ${main_prog} tag -W -T ${min_tag_len} -t ${n_thread} -O ${run} \
                -A ${linker_a} -B ${linker_b}\
                ${data_dir}/${r1_fastq} ${data_dir}/${r2_fastq} \
                2>> ${log_file}
        fi
    fi

    # Report linker detection completion
    echo -e "`date` --- Linker detection completed ----\n" >> ${log_file}


    ## Get the statistics
    # Report statistics start
    echo "
    `date` --- Statistics started  ----
    " >> ${log_file}


    if [ ${linker_b} == "none" ]; then
        # Statistics
        ${main_prog} stat -s -p -T ${min_tag_len} -t ${n_thread} ${run}.cpu \
            2>> ${log_file} 1> ${run}.stat
    else
        ${main_prog} stat -d -p -T ${min_tag_len} -t ${n_thread} ${run}.cpu \
            2>> ${log_file} 1> ${run}.stat    
    fi

    # Report statistics completion
    echo -e "`date` --- Statistics completed  ----\n" >> ${log_file}

elif [ ${experiment_type} == 'HiChIP' ]; then
    # Linker filtering for HiChIP data
    ${dep_dir}/python ${bin_dir}/util/scripts/filter_hichip_linker.py \
        --r1_file  ${data_dir}/${r1_fastq} \
        --r2_file  ${data_dir}/${r2_fastq} \
        --run  ${run} \
        --linker  ${linker_a} \
        --min_tag_len  ${min_tag_len}
    
    # Write linker filtering stats
    bash ${bin_dir}/util/scripts/write_hichip_linker_stats.sh -c ${conf}

elif [ ${experiment_type} == 'PLAC-seq' ]; then
    # Linker filtering for PLAC-seq data
    ${dep_dir}/python ${bin_dir}/util/scripts/filter_hichip_linker.py \
        --r1_file  ${data_dir}/${r1_fastq} \
        --r2_file  ${data_dir}/${r2_fastq} \
        --run  ${run} \
        --linker  ${linker_a} \
        --min_tag_len  ${min_tag_len}
    
    # Write linker filtering stats
    bash ${bin_dir}/util/scripts/write_hichip_linker_stats.sh -c ${conf}
    
fi

## Compress files
pigz -p ${n_thread} ${run}.singlelinker.paired.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.none.fastq 2>> ${log_file}

pigz -p ${n_thread} ${run}.singlelinker.single.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.conflict.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.tied.fastq 2>> ${log_file}


## Report script completion
echo -e "$0 done \n" >> ${log_file}
echo "`date`" >> ${log_file}
