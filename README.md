## Config file for ChIA PIPE

### The directory containing the executables for ChIA PIPE
    bin_dir="/projects/capurd/chia_pet/chia_pet_tool_2"

### The name of the sequencing run
    run="LHH0066"

### The type of sequencing run:
####    "miseq" - around 30 million reads
####    "hiseq" - around 300 million reads
####    "pooled" - around 1 billion reads
    run_type="miseq"

### The factor for which the IP was performed
    ip_factor="CTCF"

### Cell type
    cell_type="HG00731"

### The directory containing the input FASTQ files
    data_dir="/projects/ruan-lab/processing/fastq/${run}"

### The names of the FASTQ files
    r1_fastq="merged_${run}_R1.fastq.gz"
    r2_fastq="merged_${run}_R2.fastq.gz"


### The name of the primary genome
#### For example: "hg19", "hg38", "dm3", "mm9", "mm10"
    genome="hg38"

### The reference genome FASTA file for aligning the reads
#### (The same directory must also contain the BWA index files)
    fasta="/projects/ruan-lab/processing/genomes/hg38/hg38.fa"

### The chrom.sizes file from UCSC Genome browser
#### for the relevant genome build
    chrom_sizes="/projects/ruan-lab/processing/genomes/hg38/hg38.chrom.sizes"

### The BAM file for the ChIP-seq input control
#### (Required for spp; not required for macs2)
#### If not available, set to "none"
    input_control="none"

### The peak-calling algorithm ("macs2" or "spp")
    peak_caller="macs2"

### The folder in BASIC browser to which to upload the tracks
    basic_folder="New user testing"




## Advanced options

### The phased SNP file for allele-specific analysis
#### (If not available, set to "none")
    snp_file="none"

### How the pipeline should be executed
#### true: execute all steps in the pipeline
#### false: execute only the step submitted (mostly used for debugging)
    all_steps=true

### Should unessential intermediate files be deleted?
#### true: delete intermediate files
#### false: retain intermediate files (will take up a lot more disk space)
    clean=true

### The ChIA-PET Utilities program
    main_prog="${bin_dir}/cpu-dir/cpu-dir/cpu"

### The Juicer executable
    juicer="${bin_dir}/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar"
