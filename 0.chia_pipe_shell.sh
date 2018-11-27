#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -l mem=2GB
#PBS -l vmem=2GB
#PBS -j oe


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

# Add dependency dir to path
workdir=$( pwd)
dep_dir=$( cd ${dep_dir} && pwd )
export PATH=${dep_dir}:${PATH}
cd ${workdir}

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
    
###############
############### 1. Linker filtering

# Linker detection parameters
linker_a="none"
linker_b="none"
min_tag_len=18

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


## Compress files
pigz -p ${n_thread} ${run}.singlelinker.paired.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.none.fastq 2>> ${log_file}

pigz -p ${n_thread} ${run}.singlelinker.single.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.conflict.fastq 2>> ${log_file}
pigz -p ${n_thread} ${run}.tied.fastq 2>> ${log_file}


## Report script completion
echo -e "$0 done \n" >> ${log_file}
echo "`date`" >> ${log_file}



###############
############### 2a. Map "No linker" reads

# Setting clustering parameters
self_bp=8000
exten_bp=500

tag_name="none"
map_qual=30
suffix="UU"

## Create name of the log file
log_file=2.${run}.map_${tag_name}.log

## Print arguments to ensure correct parsing
echo "
Arguments:
    main_prog=${main_prog}
    run=${run}
    tag_name=${tag_name}
    genome=${genome}
    ctrl_genome=${ctrl_genome}
    fasta=${fasta}
    juicer=${juicer}
    out_dir=${out_dir}
" >> ${log_file}


#-- perform hybrid bwa-mem and bwa-aln mapping, 
# de-duplication, span computation, and tag clustering --#

## Perform mapping using memaln (hybrid of bwa-mem and bwa-aln)
# Report mapping start
echo "
`date` --- Mapping started for: ---
    ${run} ${tag_name}
" >> ${log_file}


${main_prog} memaln -T ${map_qual} -t ${n_thread} ${fasta} \
    ${run}.${tag_name}.fastq.gz 1> ${run}.${tag_name}.sam 2>> ${log_file}

# Compress files
pigz -p ${n_thread} ${run}.${tag_name}.sam >> ${log_file}

# Report mapping completion
echo -e "`date` --- Mapping completed ---\n" >> ${log_file}

## Pair the tags
# Report pairing start
echo -e "`date` --- Pairing paired tags ---\n" >> ${log_file}

# Pairing
${main_prog} pair -S -q ${map_qual} -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.sam.gz \
    1>${run}.${tag_name}.stat.xls 2>> ${log_file}

# Report pairing completion
echo -e "`date` --- ENDED ${run} cpu pair ---\n" >> ${log_file}

## Compute the span of the paired tags
# Report span computation start
echo -e "`date` --- Computing span of paired tags ---\n" >> ${log_file}

# Span computation
${main_prog} span -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam 2>> ${log_file} \
    1>${run}.${tag_name}.${suffix}.span.xls

# Report span computation completion
echo -e "`date` --- ENDED ${run} span pair --\n" >> ${log_file}


## Deduplicate the paired tags
# Report tag deduplication start
echo -e "`date` --- De-duplicating paired tags ${suffix} ---\n" >> ${log_file}

# Tag deduplication
${main_prog} dedup -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam \
    1> ${run}.${tag_name}.${suffix}.dedup.lc 2>> ${log_file}

# Remove intermediary file
rm ${run}.${tag_name}.${suffix}.cpu.dedup 2>> ${log_file}

# Report tag deduplication completion
echo -e "`date` --- ENDED ${run} cpu dedup ---" >> ${log_file}

## Compute the span of the non-redundant tags
# Report non-redundant span computation start
echo -e "`date` --- Computing span of paired tags ${suffix} nr ---\n" \
    >> ${log_file}

# Non-redundant span computation
${main_prog} span -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.nr.bam \
    2>> ${log_file} 1>${run}.${tag_name}.${suffix}.nr.span.xls

# Report non-redundant span computation completion
echo -e "`date` --- ENDED ${run} cpu dedup span ---\n" >> ${log_file}




###############
############### 2b. Map "Linker, one tag" reads

# Setting clustering parameters
self_bp=8000
exten_bp=500

tag_name="singlelinker.single"
map_qual=10
suffix="UxxU"

## Create name of the log file
log_file=2.${run}.map_${tag_name}.log

## Print arguments to ensure correct parsing
echo "
Arguments:
    main_prog=${main_prog}
    run=${run}
    tag_name=${tag_name}
    genome=${genome}
    ctrl_genome=${ctrl_genome}
    fasta=${fasta}
    juicer=${juicer}
    out_dir=${out_dir}
" >> ${log_file}


#-- perform hybrid bwa-mem and bwa-aln mapping, 
# de-duplication, span computation, and tag clustering --#

## Perform mapping using memaln (hybrid of bwa-mem and bwa-aln)
# Report mapping start
echo "
`date` --- Mapping started for: ---
    ${run} ${tag_name}
" >> ${log_file}


${main_prog} memaln -T ${map_qual} -t ${n_thread} ${fasta} \
    ${run}.${tag_name}.fastq.gz 1> ${run}.${tag_name}.sam 2>> ${log_file}

# Compress files
pigz -p ${n_thread} ${run}.${tag_name}.sam >> ${log_file}

# Report mapping completion
echo -e "`date` --- Mapping completed ---\n" >> ${log_file}

## Pair the tags
# Report pairing start
echo -e "`date` --- Pairing paired tags ---\n" >> ${log_file}

# Pairing
${main_prog} pair -S -q ${map_qual} -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.sam.gz \
    1>${run}.${tag_name}.stat.xls 2>> ${log_file}

# Report pairing completion
echo -e "`date` --- ENDED ${run} cpu pair ---\n" >> ${log_file}

## Compute the span of the paired tags
# Report span computation start
echo -e "`date` --- Computing span of paired tags ---\n" >> ${log_file}

# Span computation
${main_prog} span -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam 2>> ${log_file} \
    1>${run}.${tag_name}.${suffix}.span.xls

# Report span computation completion
echo -e "`date` --- ENDED ${run} span pair --\n" >> ${log_file}


## Deduplicate the paired tags
# Report tag deduplication start
echo -e "`date` --- De-duplicating paired tags ${suffix} ---\n" >> ${log_file}

# Tag deduplication
${main_prog} dedup -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam \
    1> ${run}.${tag_name}.${suffix}.dedup.lc 2>> ${log_file}

# Remove intermediary file
rm ${run}.${tag_name}.${suffix}.cpu.dedup 2>> ${log_file}

# Report tag deduplication completion
echo -e "`date` --- ENDED ${run} cpu dedup ---" >> ${log_file}

## Compute the span of the non-redundant tags
# Report non-redundant span computation start
echo -e "`date` --- Computing span of paired tags ${suffix} nr ---\n" \
    >> ${log_file}

# Non-redundant span computation
${main_prog} span -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.nr.bam \
    2>> ${log_file} 1>${run}.${tag_name}.${suffix}.nr.span.xls

# Report non-redundant span computation completion
echo -e "`date` --- ENDED ${run} cpu dedup span ---\n" >> ${log_file}



###############
############### 2c. Map "Linker, two tag" reads

# Setting clustering parameters
self_bp=8000
exten_bp=500

tag_name="singlelinker.paired"
map_qual=30
suffix="UU"

## Create name of the log file
log_file=2.${run}.map_${tag_name}.log

## Print arguments to ensure correct parsing
echo "
Arguments:
    main_prog=${main_prog}
    run=${run}
    tag_name=${tag_name}
    genome=${genome}
    ctrl_genome=${ctrl_genome}
    fasta=${fasta}
    juicer=${juicer}
    out_dir=${out_dir}
" >> ${log_file}


#-- perform hybrid bwa-mem and bwa-aln mapping, 
# de-duplication, span computation, and tag clustering --#

## Perform mapping using memaln (hybrid of bwa-mem and bwa-aln)
# Report mapping start
echo "
`date` --- Mapping started for: ---
    ${run} ${tag_name}
" >> ${log_file}


${main_prog} memaln -T ${map_qual} -t ${n_thread} ${fasta} \
    ${run}.${tag_name}.fastq.gz 1> ${run}.${tag_name}.sam 2>> ${log_file}

# Compress files
pigz -p ${n_thread} ${run}.${tag_name}.sam >> ${log_file}

# Report mapping completion
echo -e "`date` --- Mapping completed ---\n" >> ${log_file}

## Pair the tags
# Report pairing start
echo -e "`date` --- Pairing paired tags ---\n" >> ${log_file}

# Pairing
${main_prog} pair -S -q ${map_qual} -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.sam.gz \
    1>${run}.${tag_name}.stat.xls 2>> ${log_file}

# Report pairing completion
echo -e "`date` --- ENDED ${run} cpu pair ---\n" >> ${log_file}

## Compute the span of the paired tags
# Report span computation start
echo -e "`date` --- Computing span of paired tags ---\n" >> ${log_file}

# Span computation
${main_prog} span -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam 2>> ${log_file} \
    1>${run}.${tag_name}.${suffix}.span.xls

# Report span computation completion
echo -e "`date` --- ENDED ${run} span pair --\n" >> ${log_file}


## Deduplicate the paired tags
# Report tag deduplication start
echo -e "`date` --- De-duplicating paired tags ${suffix} ---\n" >> ${log_file}

# Tag deduplication
${main_prog} dedup -g -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.bam \
    1> ${run}.${tag_name}.${suffix}.dedup.lc 2>> ${log_file}

# Remove intermediary file
rm ${run}.${tag_name}.${suffix}.cpu.dedup 2>> ${log_file}

# Report tag deduplication completion
echo -e "`date` --- ENDED ${run} cpu dedup ---" >> ${log_file}

## Compute the span of the non-redundant tags
# Report non-redundant span computation start
echo -e "`date` --- Computing span of paired tags ${suffix} nr ---\n" \
    >> ${log_file}

# Non-redundant span computation
${main_prog} span -t ${n_thread} -s ${self_bp} \
    ${run}.${tag_name}.${suffix}.nr.bam \
    2>> ${log_file} 1>${run}.${tag_name}.${suffix}.nr.span.xls

# Report non-redundant span computation completion
echo -e "`date` --- ENDED ${run} cpu dedup span ---\n" >> ${log_file}


# Cluster tags using 500bp extension 
echo -e "`date` --- STARTED ${run} clustering with ${extbp} bp"\
    "extension from each side --- \n" >> ${log_file}

${main_prog} cluster -m -s ${self_bp} -B 1000 -5 5,0 -3 3,${exten_bp} \
    -t ${n_thread} -j -x -v 1 -g  -O ${run}.e500 \
    ${run}.${tag_name}.${suffix}.nr.bam 1>> ${log_file} 2>> ${log_file}

echo -e "`date` --- ENDED ${run} cpu clustering --- \n" >> ${log_file}


# Rename loop files appropriately
mv ${run}.e500.clusters.cis.chiasig.gz ${run}.e500.clusters.cis.gz
mv ${run}.e500.clusters.trans.chiasig.gz ${run}.e500.clusters.trans.gz

# Make subset file with intrachrom loops with PET_count >= 2
# for browsing in Excel
cis_file="${run}.e500.clusters.cis.gz"
be3_file="${run}.e500.clusters.cis.BE3"
zcat ${cis_file} | awk '{ if ( $7 >= 3 ) print }' > ${be3_file}


### Make .hic file for Juicebox
### BAM --> pairs --> hic
# BAM to pairs
paired_tag_bam="${run}.singlelinker.paired.UU.nr.bam"
${bin_dir}/util/pairix_src/util/bam2pairs/bam2pairs -c ${chrom_sizes} \
    ${paired_tag_bam} ${run}

# Pairs to .hic
juicer="${bin_dir}/util/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar"
hic_file="ChIA-PET_${genome}_${cell_type}_${ip_factor}_${run}_${run_type}_pairs.hic"

java -Xmx2g -jar ${juicer} pre -r \
    2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 \
    ${run}.bsorted.pairs.gz ${hic_file} ${chrom_sizes}

    
###############
############### 3. Call peaks

# Create name of the log file
log_file=3.${run}.call_peaks.log


# Print arguments to ensure correct parsing
echo "
Arguments:
    run=${run}
    peak_caller=${peak_caller}
    input_control=${input_control}
    bin_dir=${bin_dir}
    out_dir=${out_dir}
" >> ${log_file}


## Perform ChIP-Seq peak calling and plot density from ChIA-PET data

# Convert to bedgraph
# Sort bam for samtools counting bases and for BASIC visualization
if [ ! -f ${run}.for.BROWSER.bam ]
then
    echo -e "`date` Converting file formats..\n" >> ${log_file}
    samtools sort -@ 16 ${run}.singlelinker.paired.UU.nr.bam \
        -o ${run}.singlelinker.paired.UU.nr.sorted.bam
    samtools sort -@ 16 ${run}.singlelinker.single.UxxU.nr.bam \
        -o ${run}.singlelinker.single.UxxU.nr.sorted.bam
    samtools sort -@ 16 ${run}.none.UU.nr.bam \
        -o ${run}.none.UU.nr.sorted.bam

    samtools merge ${run}.for.BROWSER.bam \
        ${run}.singlelinker.paired.UU.nr.sorted.bam \
        ${run}.singlelinker.single.UxxU.nr.sorted.bam \
        ${run}.none.UU.nr.sorted.bam
    
    # Create BAM index
    samtools index ${run}.for.BROWSER.bam ${run}.for.BROWSER.bam.bai
    
    # Make bedgraph
    bedtools genomecov -ibam ${run}.for.BROWSER.bam \
        -bg > ${run}.for.BROWSER.bedgraph
    
    # Sort bedgraph
    ${bin_dir}/util/scripts/bedSort \
        ${run}.for.BROWSER.bedgraph \
        ${run}.for.BROWSER.sorted.bedgraph
    
    # Make bigwig
    ${bin_dir}/util/scripts/bedGraphToBigWig \
        ${run}.for.BROWSER.sorted.bedgraph \
        ${chrom_sizes} \
        ${run}.for.BROWSER.bigwig    
fi


if [ ${clean} == true ]
then
    ## Remove redundant BAMs
    # BAMs of reads with no linker
    rm -f ${run}.none.*.bam

    # BAMs of reads with only one usable tag
    rm -f ${run}.singlelinker.single.*.bam

    # For reads with two tags, save only the final
    # deduplicated and sorted BAM file
    ls ${run}.singlelinker.paired.*.bam | grep -v "sorted" | xargs rm -f

    # Remove FASTQ files of linker filtering
    rm -f ${run}.*.fastq.gz

    # Remove SAM files
    rm -f ${run}.*.sam.gz
    
    # Remove unsorted bedgraph
    rm -f ${run}.for.BROWSER.bedgraph
        
fi



if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
then
    # Load R
    #module load R/3.2.1
    
    # Report SPP to log file
    echo -e "`date` Peak calling using SPP..\n" >> ${log_file}
    
    # Call peaks using SPP
    z_thresh=6
    R --vanilla < ${bin_dir}/util/scripts/spp.R --args ${run}.for.BROWSER.bam \
        ${input_control} ${bin_dir} ${z_thresh}
else
    # Load MACS2
    #module load MACS/2.1.0.20151222
    
    # Report MACS2 to log file
    echo -e "`date` Peak calling using MACS2..\n" >> ${log_file}
    
    if [ ${input_control} == 'none' ] || [ ${input_control} == 'None' ]
    then
        # Call peaks using MACS2 without input control
        macs2 callpeak --keep-dup all --nomodel -t ${run}.for.BROWSER.bam \
            -f BAM -g hs -n ${run}.no_input_all 1>> ${log_file} 2>> ${log_file}
    else
        # Call peaks using MACS2 with input control
        macs2 callpeak --keep-dup all --nomodel -t ${run}.for.BROWSER.bam \
            -c ${input_control} \
            -f BAM -g hs -n ${run}.all 1>> ${log_file} 2>> ${log_file}
    fi
fi

echo -e "`date` ENDED ${run} peak calling ..\n" >> ${log_file}
echo "$0 done" >> ${log_file}


#### Annotate loops with peak support and call CCDs
# Annotate loops
be3_file="${run}.e500.clusters.cis.BE3"
peak_file="${run}.for.BROWSER.spp.z6.broadPeak"

${bin_dir}/util/scripts/annotate_loops_with_peak_support.py \
    -l ${be3_file} -p ${peak_file}

# Create subsetted files by peak support
# 2 anchors
cat ${be3_file}.peak_annot | awk '{ if ( $8 == 2 ) print }' \
    > ${be3_file}.peak_annot.E2

# 1 or more anchor
cat ${be3_file}.peak_annot | awk '{ if ( $8 >= 1 ) print }' \
    > ${be3_file}.peak_annot.BE1

####
# Sort peak-supported loops by PET count
#sort -k7,7n ${be3_file}.peak_annot.E2 > ${be3_file}.peak_annot.E2.freq_sorted

## Get PET count cutoff for calling CCDs
## Top one-third (67th percentile)

# Line number in sorted file
cutoff_line=$( echo -e \
    "$( cat ${be3_file}.peak_annot.E2.freq_sorted | wc -l ) * 0.67 / 1" | bc )

# PET count cutoff
min_pet_count=$( sed "${cutoff_line}q;d" \
    ${be3_file}.peak_annot.E2.freq_sorted | awk '{ print $7 }' )

### Call CCDs
#module load bedtools/2.26.0

# Set input file for calling CCDs
ccd_input="${be3_file}.peak_annot.E2"

# Make bed file of loop spans
#awk -v OFS='\t' '{ print $1, $3, $5, $7 }' ${ccd_input} \
#    > ${ccd_input}.inner

${bin_dir}/util/scripts/get_loop_anchor_midpoints.py -l ${ccd_input} \
    -m ${min_pet_count}


# Sort loops
sort -k1,1 -k2,2n ${ccd_input}.anchor_mids > ${ccd_input}.anchor_mids.sorted

# Merge loops
bedtools merge -i ${ccd_input}.anchor_mids.sorted \
    > ${ccd_input}.anchor_mids.sorted.merge

# Add ccd span
awk -v OFS='\t' '{ $4 = $3 - $2} 1' ${ccd_input}.anchor_mids.sorted.merge \
    > ${ccd_input}.anchor_mids.sorted.merge.span

# Filter out CCDs with span < 10kb
awk -v OFS='\t' '{ if ( $4 > 25000 ) print $0 }' \
    ${ccd_input}.anchor_mids.sorted.merge.span \
    > ${ccd_input}.anchor_mids.sorted.merge.span.ccd

# Remove CCD intermediate files
rm -f *.anchor_mids.sorted 
rm -f *.anchor_mids.sorted.merge
rm -f *.anchor_mids.sorted.merge.span



###############
############### 4. Extract summary stats

## Set the output file
out_file=${run}.final_stats.tsv
rm -f ${out_file}

## Get library ID
echo -e "Library_ID\t"${run} >> ${out_file}

# Get library type
echo -e "Library_type\t"${run_type} >> ${out_file}

# Get reference genome
echo -e "Reference_genome\t"${genome} >> ${out_file}

# Get cell type
echo -e "Cell_type\t"${cell_type} >> ${out_file}

# Get IP factor
echo -e "Factor\t"${ip_factor} >> ${out_file}

# Create contact-map URL
hic_file="ChIA-PET_${genome}_${cell_type}_${ip_factor}_${run}_${run_type}_pairs.hic"
url="http://ctencode01.jax.org/chiapet/dev/${hic_file}"

echo -e "Contact-map_URL\t"${url} >> ${out_file}

## PET count
# Get PET count
n_read_pair=$( cat ${run}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )

## Get linker statistics
read_pair_link=$( cat ${run}.stat | grep "Linker detected" | \
    awk -F '[ \t]' '{print $3}' )

frac_link=$( echo -e "${read_pair_link} / ${n_read_pair}" | bc -l | xargs printf "%.2f\n")

# Write PET count
n_read_pair=$( printf "%'.f\n" ${n_read_pair} )
echo -e "Total_read_pairs\t"${n_read_pair} >> ${out_file}

# Write linker statistics
read_pair_link=$( printf "%'.f\n" ${read_pair_link} )
echo -e "Read_pairs_with_linker\t"${read_pair_link} >> ${out_file}
echo -e "Fraction_read_pairs_with_linker\t"${frac_link} >> ${out_file}

# Write one tag vs two tag
one_tag=$( grep "Single Linker 1 tag (SL/ls)" ${run}.stat | cut -f2 )
two_tag=$( grep "Single Linker 2 tags (SL/ls)" ${run}.stat | cut -f2 )

one_tag=$( printf "%'.f\n" ${one_tag} )
two_tag=$( printf "%'.f\n" ${two_tag} )

echo -e "One_tag\t"${one_tag} >> ${out_file}
echo -e "PET\t"${two_tag} >> ${out_file}

## Mapping
# Get uniquely mapped PET count 
unique=$( cat ${run}.singlelinker.paired.UU.span.xls | grep "Total pairs" | \
    awk -F '[\t]' '{print $2}' )

# Get uniquely mapped and non-redundant PET count 
nr=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | grep "Total pairs" | \
    awk -F '[\t]' '{print $2}' )

# Compute redundancy
redun=$( echo "(${unique} - ${nr}) / ${unique}" | bc -l )

# Write uniquely mapped PET count
unique=$( printf "%'.f" ${unique} )
echo -e "Uniquely_mapped_PET\t"${unique} >> ${out_file}

# Write unique mapped and non-redundant PET count
nr=$( printf "%'.f" ${nr} )
echo -e "Non-redundant_PET\t"${nr} >> ${out_file}

# Write redundancy
redun=$( printf %.2f ${redun} )
echo -e "Redundancy\t"${redun} >> ${out_file}

# Write non-redundant tags
nr_tag=$( samtools view -c ${run}.for.BROWSER.bam )
echo -e "Non-redundant_tag\t"${nr_tag} >> ${out_file}

## Get number of peaks
if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
then
    n_peak=$( cat ${run}.for.BROWSER.spp.z6.broadPeak | wc -l )
else
    if [ ${input_control} == "none" ]
    then
        n_peak=$( cat ${run}.no_input_all_peaks.narrowPeak | wc -l )
    else
        n_peak=$( cat ${run}.all_peaks.narrowPeak | wc -l )
    fi
fi

n_peak=$( printf "%'.f" ${n_peak} )
echo -e "Peak\t"$n_peak >> ${out_file}


## Interaction types
# Get self-ligation PET count
self_lig=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==4)print $2}' )

self_lig=$( printf "%'.f" ${self_lig} )
echo -e "Self-ligation_PET\t"${self_lig} >> ${out_file}

# Get inter-ligation PET count (intra-chr)
intra_chr_pet=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==5)print $2}' )


# Get inter-ligation PET count (inter-chr)
inter_chr_pet=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==2)print $2}' )

# Compute ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( echo "${intra_chr_pet} / ${inter_chr_pet}" | bc -l )

# Compute inter-ligation PET count (all)
inter_lig_all=$( echo "${intra_chr_pet} + ${inter_chr_pet}" | bc )

# Write inter-ligation PET count (all)
inter_lig_all=$( printf "%'.f" ${inter_lig_all} )
echo -e "Inter-ligation_PET\t"${inter_lig_all} >> ${out_file}

# Write inter-ligation PET count (intra-chr)
intra_chr_pet=$( printf "%'.f" ${intra_chr_pet} )
echo -e "Intra-chr_PET\t"${intra_chr_pet} >> ${out_file}

# Write inter-ligation PET count (inter-chr)
inter_chr_pet=$( printf "%'.f" ${inter_chr_pet} )
echo -e "Inter-chr_PET\t"${inter_chr_pet} >> ${out_file}

# Write ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( printf %.2f ${pet_ratio} )
echo -e "ratio_of_intra/inter_PET\t"${pet_ratio} >> ${out_file}


## Singleton
# Get singleton PET count (all)
singleton=$(zcat *clusters*.gz | awk '$7==1{print}' | wc -l)
singleton=$( printf "%'.f" ${singleton} )
echo -e "Singleton\t"$singleton >> ${out_file}

# Get singleton PET count (intra-chr)
intra_singleton=$(zcat *cis.gz | awk '$7==1{print}' | wc -l)
intra_singleton=$( printf "%'.f" ${intra_singleton} )
echo -e "Intra-chr_singleton\t"$intra_singleton >> ${out_file}

# Get singleton PET count (inter-chr)
inter_singleton=$(zcat *trans.gz | awk '$7==1{print}' | wc -l)
inter_singleton=$( printf "%'.f" ${inter_singleton} )
echo -e "Inter-chr_singleton\t"$inter_singleton >> ${out_file}


## Clusters (overall)
# Get cluster count
total_cluster_number=$(zcat *clusters*.gz | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}

# Get intra-chr cluster count
intra_cluster=$( zcat *cis.gz | awk '$7 >=2 {print}' | wc -l )

# Get inter-chr cluster count
inter_cluster=$( zcat *trans.gz | awk '$7 >=2 {print}' | wc -l)


# Compute ratio of intra-chr to inter-chr clusters
cluster_ratio=$( echo "${intra_cluster} / ${inter_cluster}" | bc -l )
cluster_ratio=$( printf %.2f ${cluster_ratio} )

# Write cluster ratio
echo -e "ratio_of_intra/inter_cluster\t"${cluster_ratio} >> ${out_file}


## Clusters (intra-chr)

# Write intra-chr cluster count
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}

# Get intra-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
	intra_pets_number=$(zcat *cis.gz | \
	    awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
	    xargs printf "%'.f")
	
	echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done

# Get intra-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *cis.gz | awk '$7 >10 {print}' | \
    wc -l | xargs printf "%'.f") >> ${out_file}

## Write loops with peak support
#if [ ${ip_factor} == 'CTCF' ]; then
#    loops_pk_supp=$( cat ${run}.e500.clusters.cis.BE3.peak_annot.E2 | wc -l )
#else
#    loops_pk_supp=$( cat ${run}.e500.clusters.cis.BE3.peak_annot.E2 | wc -l )
#fi
#
#echo -e "Loops_with_peak_support\t"${loops_pk_supp} >> ${out_file}

## Clusters (inter-chr)
# Write inter-chr cluster count
inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}


# Get inter-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
	inter_pets_number=$(zcat *trans.gz | \
	    awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
	    xargs printf "%'.f")
	
	echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done

# Get inter-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *trans.gz | \
    awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}


###############
############### 5. Phase loops

if [ ${all_steps} == true ] && [ ${snp_file} != "none" ]
then
    ## Create name of the log file
    log_file=5a.${run}.phase_loops.log

    # Print arguments to ensure correct parsing
    echo "
    Arguments:
        run=${run}
        fasta=${fasta}
        snp_file=${snp_file}
        out_dir=${out_dir}
        bin_dir=${bin_dir}
    " >> ${log_file}


    # Load the required modules
    #module load R/3.2.1
    #module load perl/5.24.0
    #module load bedtools/2.26.0

    # Create file names
    bam_file=${run}.for.BROWSER.bam
    loop_file=${run}.e500.clusters.cis.BE3
    pileup_file=${run}.total_snp_pileup.txt

    if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
    then
        peak_file=${run}.for.BROWSER.spp.z6.broadPeak
    else
        peak_file=${run}.no_input_all_peaks.narrowPeak
    fi


    # Separate the SNP file into one file for each chromosome
    # (named after the chromosome)
    echo -e "`date` Splitting SNP file by chromosome..\n" >> ${log_file}

    awk '{ print >> "temp_"$1 }' ${snp_file}


    # Count the SNP alleles for each chromosome
    echo -e "`date` Counting SNP alleles on each chromosome..\n" >> ${log_file}

    ### Get chromosome names from the fasta file
    chroms=$(grep ">" ${fasta} | sed -e "s/>//g")
    echo -e "`date` Detected the following chroms in the FASTA file..\n" \
        >> ${log_file}

    for chrom in ${chroms}
    do
        echo -e "\t${chrom}" >> ${log_file}
    done


    ### Pileup chromosomes
    echo -e "`date` Performing mpileup on chromosome:\n" \
        >> ${log_file}

    job_ids=()

    for chrom in ${chroms}
    do
        # Skip chrM
        if [ ${chrom} == 'chrM' ] || [ ${chrom} == 'chrY' ]
        then
            echo -e "\t (Skipping ${chrom})" >> ${log_file}
            continue
        fi
    
        # Create a separate PBS job for each chromosome
        # to count SNP alleles with samtools
        echo -e "\t ${chrom}" >> ${log_file}
        job_name="temp_job_${chrom}"
    
        # Change to output directory and load modules
        pbs="cd ${out_dir};"
        pbs+="module load samtools/0.1.19;"
    
        # Samtools command
        # By default, samtools/0.1.19 has BAQ filtering ON
        # It is important to retain this if changing to a different version
        pbs+="samtools mpileup -f ${fasta} "
        pbs+="-l temp_${chrom} ${bam_file} "
        pbs+="> temp_counts.${chrom}.mileup"
    
        # Submit job
        jid=$( echo -e ${pbs} 2>> ${log_file} | \
            qsub -l nodes=1:ppn=20,mem=60gb,walltime=8:00:00 -N ${job_name})
        sleep 0.5
        job_ids+=" ${jid}"
    done


    # Wait for all chromosome jobs to finish
    echo -e "`date` Waiting for job completion... \n" >> ${log_file}

    for jid in ${job_ids}
    do
        echo -e "\t${jid}" >> ${log_file}
    
        while [ -n "$( qstat | grep ${jid} )" ]
        do
            sleep 1
        done
    done

    # Report completion of counting SNP alleles
    echo -e "`date` Completed counting SNP alleles..\n" >> ${log_file}


    ### Merge SNP allele counts from all chromosomes
    echo -e "`date` Merging SNP allele counts from all chroms..\n" \
        >> ${log_file}

    cat *mileup | awk 'BEGIN{FS="\t";OFS="\t"}{chrom=$1;pos=$2;ref=$3;\
    num=$4+$7+$10;bases=$5""$8""$11;qual=$6""$9""$12;\
    print chrom,pos,ref,num,bases,qual}' | awk 'NF==6{print}' \
        > ${pileup_file}

    echo -e "`date` Completed merging SNP allele counts..\n" >> ${log_file}


    # Remove temporary files
    echo -e "`date` Removing temp files..\n" >> ${log_file}

    #rm temp_*

    echo -e "`date` Completed removing temp files..\n" >> ${log_file}

    ### Count phased SNPs
    # Do binomial test to identify significantly biased SNPs
    # Output file name: phased_snp_qvalue.txt
    echo -e "`date` Counting SNP allele bias and assessing significance.." \ 
        >> ${log_file}

    # Allele-biased SNPs (R)
    Rscript ${bin_dir}/compute_phased_SNP.R \
        ${run} ${pileup_file} ${snp_file} TRUE 0.1

    # Create BED file of All SNPs
    cat ${run}.snp_qvalues.txt | grep -v biased | awk -v OFS='\t' \
        '{ print $1, $2, $2+200, $7, "+", $4, $5, $6, "153,153,153"}' \
        > ${run}.snp_qvalues.bed

    # Create BED file of Phased SNPs
    cat ${run}.snp_qvalues.txt | grep Yes | awk -v OFS='\t' \
        '{ 
        if ($4 >= $5) 
            print $1, $2, $2+200, $7, "+", $4, $5, $6, "26,152,80"
        else
            print $1, $2, $2+200, $7, "+", $4, $5, $6, "215,48,39"
        }' \
        > ${run}.snp_qvalues.phased.bed


    # Make bedgraph of Paternal allele counts
    cat phased_snp_qvalue.txt | grep -v biased | \
        awk -v OFS='\t' '{ print $1, $2, $2+100, $4 }' \
        > ${run}.snp_coverage.Paternal.100.bedgraph

    # Make bedgraph of maternal allele counts
    cat phased_snp_qvalue.txt | grep -v biased | \
        awk -v OFS='\t' '{ print $1, $2, $2+10, $5 }' \
        > ${run}.snp_coverage.Maternal.10.bedgraph


    echo -e "`date` Completed counting SNP allele bias..\n" >> ${log_file}


    #######
    cat phased_snp_qvalue.Yes.txt | grep -v biased | \
        awk -v OFS='\t' '{ print $1, $2, $2+1, $4 }' \
        > ${run}.snp_coverage.Yes.Paternal.1.bedgraph
    ########


    ### Compute phased interactions
    # phased_snp_qvalue.txx is generated by last step
    # FIXME: unique_nodes_interaction could be output file from cpu
    # (non-overlapping)
    echo "`date` Identifying biased loops (from raw loops and biased SNPs)..
    " >> ${log_file}


    perl ${bin_dir}/compute_phased_loops.pl \
        phased_snp_qvalue.txt ${peak_file} ${loop_file}


    # Create two subset BED files
    # BED file of only phased peaks
    cat ${peak_file}.SNP_Phased.browser.bed | grep -v Unphased \
        > ${peak_file}.SNP_Phased.Paternal_or_Maternal.bed

    # BED file of only unphased peaks
    cat ${peak_file}.SNP_Phased.browser.bed | grep Unphased \
        > ${peak_file}.SNP_Phased.Unphased.bed



    echo -e "`date` Completed identifying biased loops..\n" >> ${log_file}
    echo "`date` $0 done" >> ${log_file}
fi

