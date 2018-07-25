#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python

import sys
import os
import argparse
import gzip
import numpy as np


def get_column_index_of_sample_id(entry, sample_id):
    """
    Get the column index of a sample ID in the 1000 Genomes Project
    raw VCF data file.
    
    Args:
        entry (list):
            A list of a row of data from the VCF file.
        sample_id (str):
            The ID of the relevant sample.
    
    Returns:
        col_idx (int):
            The column index of the sample ID of interest in the
            VCF file.
    """                    
    # Np array for easy searching
    entry = np.array(entry)
    
    # Look up column index of sample
    try:
        col_idx = np.where(entry == sample_id)[0][0]
    except IndexError:
        # If the sample ID is not found
        print ('Warning: sample ID %s was not found in the '
               'VCF file.' % sample_id)
        sys.exit()
    
    # Print column index
    print 'Sample ID %s has column index %s' % (
        sample_id, col_idx)
    
    return col_idx


def vcf_to_phased_snp(vcf_dir, child_id, paternal_id, maternal_id, out_dir):
    """
    Convert the raw VCF data file of genotypes from the 1000 Genomes Project
    to the four-column phased SNP file needed as input to the haplotype
    section of the ChIA-PET pipeline.
    
    Args:
        vcf_dir (str):
            The directory containing the VCF files of phased genotypes from
            the 1000 Genomes Project.
        sample_id (str):
            The sample ID of the child for which to extract the data.
        paternal_id (str):
            The sample ID of the father of that child.
        sample_id (str):
            The sample ID of the mother of that child.
        out_dir (str):
            The output directory to which to write the resulting
            SNP file.
    """
    # Set the human chromosome labels
    chroms = []
    for i in range(1, 23) + ['X']:
        chroms.append('chr%s' % i)
    
    #
    ## Subset chroms when debugging
    chroms = ['chr22']
    #
    
    # Set the template of the VCF file names
    # to support iteration by chromosome
    vcf_template = ('ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.'
                    '20130502.genotypes.vcf.gz')
    vcf_template = os.path.join(vcf_dir, vcf_template)
    
    # Valid variant alleles to be retained
    # from the VCF file
    nucs = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']
    
    # Initialize output result array
    snps = []
    
    for chrom in chroms:
        with gzip.open(vcf_template % chrom, 'r') as f:
            for line in f:
                # Skip metadata lines at the top of the file
                # (indicated with double hashtag)
                if '##' in line:
                    continue
            
                entry = line.strip().split()
            
                # The primary header line has a single hash tag
                if '#' in entry[0]:
                    # Get the column index of the father
                    pat_col_idx = get_column_index_of_sample_id(
                        entry, paternal_id)
                    
                    # Get the column index of the mother
                    mat_col_idx = get_column_index_of_sample_id(
                        entry, maternal_id)
                    
                    continue
                
                # Parse the data lines
                chrom = entry[0]
                
                end = int(entry[1])
                start = str(end - 1)
                end = str(end)
                
                variant_name = entry[2]
                ref = entry[3]
                alt = entry[4]
                pat_binary_geno_str = entry[pat_col_idx].split(':')[0]
                mat_binary_geno_str = entry[pat_col_idx].split(':')[0]
                
                ## IMPORTANT
                # Only consider simple SNPs
                # Filter our complex variants
                if len(ref) != 1 or ref not in nucs:
                    continue
                
                if len(alt) !=1 or alt not in nucs:
                    continue
                
                # Standardize allele labels
                ref = ref.upper()
                alt = alt.upper()
                
                
                ## IMPORTANT
                # Only consider unambiguous SNPs
                # That is, SNPs where one parent is homozygous for the
                # reference allele and the other parent is homozygous for
                # the alternative allele.
                
                if pat_binary_geno_str != mat_binary_geno_str:
                        print pat_binary_geno_str, mat_binary_geno_str
                
                # Skip SNPs where parents are heterozygous
                if pat_binary_geno_str in ['0|1', '1|0']:
                    continue
                if mat_binary_geno_str in ['0|1', '1|0']:
                    continue
                

                # Skip SNPs where parents have the same allele
                if pat_binary_geno_str == mat_binary_geno_str:
                    continue
                
                # Get the paternal allele
                if pat_binary_geno_str == '0|0':
                    pat_allele = ref
                elif pat_binary_geno_str == '1|1':
                    pat_allele = alt
                else:
                    print ('Atypical genotype for paternal ID of variant '
                           '%s.' % variant_name)
                    continue
                
                # Get the maternal allele
                if mat_binary_geno_str == '0|0':
                    mat_allele = ref
                elif mat_binary_geno_str == '1|1':
                    mat_allele = alt
                else:
                    print ('Atypical genotype for maternal ID of variant '
                           '%s.' % variant_name)
                    continue
                
                # Make the final allele str
                child_geno_str = '|'.join([pat_allele, mat_allele])
                
                # Append results
                snps.append([chrom, start, end, child_geno_str])
    
    
    # Set output file name
    out_file = '%s.phased_snp.bed' % child_id
    out_file = os.path.join(out_dir, out_file)
    
    # Write results
    with open(out_file, 'w') as f:
        for snp in snps:
            f.write('\t'.join(snp) + '\n')            
            

def parse_command_line_args():
    """
    Parse command-line arguments.
    
    Returns:
        args (class 'argparse.Namespace'):
            An object containing the parsed command-line arguments.
            For every command-line option, the values are stored as follows:
            args.{option} 
    """
    # Initiate the argument parser
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    
    # Indicate the required arguments
    required.add_argument(
        '-v', '--vcf_dir', required=True,
        help=('The directory containing the VCF files of phased genotypes'
              'from the 1000 Genomes Project.'))

    required.add_argument(
        '-c', '--child_id', required=True,
        help=('The sample ID of the child for which to extract the data. '
              'Must be part of a 1000G trio.'))
    
    required.add_argument(
        '-p', '--paternal_id', required=True,
        help=('The sample ID of the father of that child.'))
    
    required.add_argument(
        '-m', '--maternal_id', required=True,
        help=('The sample ID of the mother of that child.'))
    
    required.add_argument(
        '-o', '--out_dir', required=True,
        help=('The output directory to which to write the resulting'
              'SNP file.'))
        
    # Parse the arguments from the command-line input
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    
    # Read command-line arguments
    args = parse_command_line_args()
    
    # Compute phased SNPs
    vcf_to_phased_snp(
        args.vcf_dir,
        args.child_id,
        args.paternal_id,
        args.maternal_id,
        args.out_dir)
