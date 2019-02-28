#!/bin/bash

# ChIA-PIPE
#         Step 2: Map the reads
# 2018
# The Jackson Laboratory for Genomic Medicine


## The help message:
function usage
{
    echo -e "usage: bash 2.map.sh --conf ${conf} --tag_name ${tag_name} 
    " 
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -t | --tag_name )       shift
                                tag_name=$1   
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

# Set mapping parameters and file suffix
if [ ${tag_name} == "singlelinker.single" ]; then
    map_qual=10
    suffix="UxxU"
else
    map_qual=30
    suffix="UU"
fi

# Update mapping parameters for hybrid data
if [ ${hybrid} != "none" ]; then
    map_qual=0
fi

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

# Mapping
${main_prog} memaln -T ${map_qual} -t ${n_thread} ${fasta} \
    ${run}.${tag_name}.fastq.gz 1> ${run}.${tag_name}.sam 2>> ${log_file}

# Compress files
pigz -p 4 ${run}.${tag_name}.sam >> ${log_file}

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


## Cluster tags and make files for visualization with Juicebox and HiGlass
if [ ${tag_name} == "singlelinker.paired" ]
then
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
    
    
    # Copy .hic file to URL directory
    # for automatic viewing of 2D contact map
    url_dir="/projects/encode/to_ctencode01/dev"
    
    if [ -d ${url_dir} ]; then
        cp ${hic_file} ${url_dir}
    fi
    
    ### Multi-resolution matrix file for HiGlass
    ## Minimum resolution for HiGlass
    #resolution=1000
    
    ## Divide chromosomes into bins
    #${bin_dir}/conda/bin/cooler makebins \
    #    ${chrom_sizes} \
    #    ${resolution} \
    #    > ${out_dir}/temp_chrom_bins.bed
    
    ## Make .cool file (1000bp resolution)
    #${bin_dir}/conda/bin/cooler cload pairix -p 20\
    #    ${out_dir}/tmp_chrom_bins.bed \
    #    ${out_dir}/${run}.e500.cooler.bsorted.pairs.gz \
    #    ${out_dir}/${run}.e500.higlass.cool
    
    ## Make final, multi-resolution .cool file (with zooming ability)
    #${bin_dir}/conda/bin/cooler zoomify -p 20\
    #    --no-balance ${out_dir}/${run}.e500.higlass.cool
    
    ## Remove intermediary file
    #rm ${out_dir}/temp_chrom_bins.bed
    #rm ${out_dir}/${run}.e500.higlass.cool

fi

# Report completion of mapping
echo "$0 done" >> 2.${run}.map_${tag_name}.done
echo "$0 done" >> ${log_file}
echo "`date`" >> ${log_file}
