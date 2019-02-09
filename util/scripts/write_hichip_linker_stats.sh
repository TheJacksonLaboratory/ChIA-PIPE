## The help message:
function usage
{
    echo -e "usage: bash write_hichip_linker_stats.sh -c CONF
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


source ${conf}

# Get total number of read pairs
pairs=$( zcat ${data_dir}/${r1_fastq} | wc -l )
pairs=$( echo "${pairs} / 4" | bc )

# Get number of read pairs with no linker
none=$( cat ${run}.none.fastq | wc -l )
none=$( echo "${none} / 8" | bc )

# Get number of read pairs with linker
linker=$( echo "${pairs} - ${none}" | bc )

# Get number of linker-1_tag read pairs
one_tag=$( cat ${run}.singlelinker.single.fastq | wc -l )
one_tag=$( echo "${one_tag} / 8" | bc )

# Get number of linker-2_tag read pairs 
two_tag=$( cat ${run}.singlelinker.paired.fastq | wc -l )
two_tag=$( echo "${two_tag} / 8" | bc )

# Write results
echo -e "Total pairs\t${pairs}" >> ${run}.stat
echo -e "Linker detected\t${linker}" >> ${run}.stat
echo -e "Single Linker 2 tags (SL/ls)\t${two_tag}" >> ${run}.stat
echo -e "Single Linker 1 tag (SL/ls)\t${one_tag}" >> ${run}.stat
