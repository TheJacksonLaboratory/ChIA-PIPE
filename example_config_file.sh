##### Config file for ChIA PIPE ####

### 1) Library information

# The name of the sequencing run
# Unique identifier for the library being processed
run="LDK0004-ds"

# Experiment type
#   "ChIA-PET"
#   "HiChIP"
#   "PLAC-seq"
experiment_type="ChIA-PET"

# The type of sequencing run:
#    "miseq" - around 30 million reads
#    "hiseq" - around 300 million reads
#    "pooled" - around 1 billion reads
run_type="miseq"

# The factor for which the IP was performed
ip_factor="RNAPII"

# Cell type
cell_type="Kc167"

# The directory containing the input FASTQ files
data_dir="../fastq/"

# The names of the FASTQ files
r1_fastq="LDK0004-ds_*_R1_*.fastq.gz"
r2_fastq="LDK0004-ds_*_R2_*.fastq.gz"



### 2) Dependency information

# The directory containing the executables for ChIA-PIPE
bin_dir="/projects/encode/chia_pipe"

# The directory of the local install of ChIA-PIPE dependencies
dep_dir="dep_dir"

# The name of the primary genome
# For example: "hg19", "hg38", "dm3", "mm9", "mm10"
genome="dm3"

# The reference genome FASTA file for aligning the reads
# (The same directory must also contain the BWA index files)
fasta="/projects/ruan-lab/processing/genomes/dm3/dm3.fa"

# The chrom.sizes file from UCSC Genome browser
# for the relevant genome build
chrom_sizes="/projects/ruan-lab/processing/genomes/dm3/dm3.chrom.sizes"

# The BAM file for the ChIP-seq input control
# (Required for spp; not required for macs2)
# If not available, set to "none"
input_control="none"

# The peak-calling algorithm ("macs2" or "spp")
peak_caller="macs2"

# The folder in BASIC browser to which to upload the tracks
basic_folder="New user testing"

# BED file of promoters for annotating loops
promoter_bed_file="none"

# BED file of enhancers for annotating loops 
enhancer_bed_file="none"


## 3) Advanced options: how to run the pipeline

# The phased SNP file for allele-specific analysis
# (If not available, set to "none")
snp_file="none"

# How the pipeline should be executed
# true: execute all steps in the pipeline
# false: execute only the step submitted (mostly used for debugging)
all_steps=true

# Should unessential intermediate files be deleted?
# true: delete intermediate files
# false: retain intermediate files (will take up a lot more disk space)
clean=true

# The ChIA-PET Utilities program
main_prog="${bin_dir}/util/cpu-dir/cpu-dir/cpu"

# The Juicer executable
juicer="${bin_dir}/util/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar"

# The number of threads to use on the HPC
n_thread=20

# The amount of memory to use on the HPC (in GB)
mem=60



## 4) Expert options: detailed parameters

# The bridge linker sequence used during ChIA-PET
linker_a="ACGCGATATCTTATCTGACT"

# The second half linker sequence (if two linker sequences were used)
linker_b="none"

# The minimum tag length required for aligning to the reference
# genome after linker filtering
min_tag_len=18

# The span threshold for calling an interaction as 
#self-ligation vs inter-ligation
self_bp=8000

# The tag extension length during loop clustering
exten_bp=500

# The Z-score threshold for peak calling with SPP
z_thresh=6

# The name of the genetic cross if working with a hybrid strain
hybrid="none"



### 5) Uncomment below to module load dependencies
# module load pigz
# module load java/1.7.0
# module load perl/5.26.0
# module load bedtools/2.26.0
# module load samtools/1.5
# module load R/3.2.1
# module load MACS/2.1.0.20151222
###
