### Config file for ChIA-PET Tool

# The name of the sequencing run
run="LHK0002"

# The type of sequencing run:
#    "miseq" - around 30 million reads
#    "nextseq" - around 300 million reads
#    "hiseq" - around 300 million reads
#    "pooled" - around 1 billion reads
run_type="miseq"

# The factor for which the IP was performed
ip_factor="CTCF"

# Cell type
cell_type="K562"

# The directory containing the input FASTQ files
data_dir="/projects/capurd/chia_pet/ENCODE_SOURCE_DATA/HUMAN/LHK0002_LKK0002N"
r1_fastq="LHK0002_S2_L001_R1_001.fastq.gz"
r2_fastq="LHK0002_S2_L001_R2_001.fastq.gz"

# The directory containing the executables for ChIA-PET Tool 
bin_dir="/projects/capurd/chia_pet/chia_pet_tool_2"

# The name of the primary genome
# "Must be one of: hg18, hg19, hg38, dMel, mm9, mm10, anasPlat1, 
#  bTaurus3, canFam3, equCab2, galGal4, Pf3D7, sacCer3, 
#  sCerS288c, susScr3, or TAIR10"
#  Otherwise, it should be the path of the chrom.sizes file
# for the relevante genome.
genome="hg38"

# The name of the spike-in control genome (if applicable)
ctrl_genome="none"

# The reference genome FASTA file for aligning the reads
# (The same directory must also contain the BWA index files)
fasta="/projects/capurd/chia_pet/genomes/hg38/hg38.fa"

# The chrom.sizes file from UCSC Genome browser
# for the relevant genome build
chrom_sizes="/projects/capurd/chia_pet/genomes/hg38/hg38.chrom.sizes"

# The Juicer executable
juicer="/projects/tjongh/chiapet_script/encode_singleLinker/\
juicer_tools.1.6.2_linux_jcuda.0.8.jar"

# The peak-calling algorithm ("macs2" or "spp")
peak_caller="macs2"

# The BAM file for the ChIP-seq input control
# (Required for spp; not required for macs2)
# If not available, set to "none"
input_control="none"

# The phased SNP file for allele-specific analysis
# (If not available, set to "none")
snp_file="none"

# How the pipeline should be executed
# true: execute all steps in the pipeline
# false: execute only the step submitted (mostly used for debugging)
all_steps=true
