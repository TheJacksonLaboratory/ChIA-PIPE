#!/bin/bash

# ChIA-PIPE
#         Step 4: Extract summary statistics from output files
# 2018
# The Jackson Laboratory for Genomic Medicine


## The help message:
function usage
{
    echo -e "usage: bash 4.extract_summary_stats.sh --conf ${conf} 
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


## Source config file
source ${conf}

# Add dependency dir to path
export PATH=${dep_dir}:${PATH}

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
