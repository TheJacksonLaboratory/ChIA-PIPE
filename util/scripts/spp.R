# Parse command-line arguments
args = commandArgs(trailingOnly=T)
chip_bam_file = args[1]
input_bam_file = args[2]
bin_dir = args[3]
z_thresh = as.numeric(args[4])

# Load packages and set library path environment variable
library("Rsamtools", lib.loc=bin_dir)
library("fastcluster", lib.loc=bin_dir)
library("spp", lib.loc=bin_dir)
.libPaths(bin_dir)


# Read the BAM files for the ChIP sample and the input-control sample
chip_reads = read.bam.tags(chip_bam_file)
input_reads = read.bam.tags(input_bam_file)


# Estimate of the binding-peak separation distance, 
# cross-correlation profile, and tag-quality bin-acceptance information
binding_char = get.binding.characteristics(
    chip_reads, srange=c(270, 1000), bin=5, accept.all.tags=T)


# Select tags based on the binding characteristics
chip_data = select.informative.tags(chip_reads, binding_char)
input_data = select.informative.tags(input_reads, binding_char)


# Remove or restrict singular positions with extremely high tag count
# relative to their neighborhood
chip_filt = remove.local.tag.anomalies(chip_data)
input_filt = remove.local.tag.anomalies(input_data)

# Get window size
window_half_size = binding_char$whs

# Identify broad regions of enrichment for the specified scale 
broad_clusters = get.broad.enrichment.clusters(
    chip_filt, input_filt, window.size=window_half_size, z.thr=z_thresh,
    tag.shift=round(binding_char$peak$x / 2))

# Write output file of broad clusters
suffix = paste0(".spp.z", z_thresh, ".broadPeak")
out_file = gsub(".bam", suffix, chip_bam_file)
write.broadpeak.info(broad_clusters, out_file)
