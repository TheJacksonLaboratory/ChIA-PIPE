#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python

import argparse
import pandas as pd
from scipy.stats import binom_test
from statsmodels.sandbox.stats.multicomp import multipletests


def read_pileup_file(pileup_file, min_read_count=1):
    """
    Read the pileup file from samtools (with all chroms combined).
    
    Args:
        pileup_file (str):
            The pileup file.
        min_read_count (int):
            The minimum number of reads a SNP must be covered by in order
            to be retained in the data frame.
    
    Returns:
        pileup (pd.DataFrame):
            The pileup data.
    """
    pileup = []
    with open(pileup_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            chrom = entry[0]
            loc = int(entry[1])
            ref = entry[2]
            num = int(entry[3])
            bases = entry[4]
            qual = entry[5]
            
            pileup.append([chrom, loc, ref, num, bases, qual])
    
    # Read pileup file
    pileup = pd.DataFrame(
        pileup,
        columns=['chrom', 'loc', 'ref', 'num', 'bases', 'qual'])
    
    # Only retain SNPs covered by at least one read
    pileup = pileup[pileup['num'] >= min_read_count]
    
    return pileup


def read_snp_file(snp_file):
    """
    Read the file of phased SNPs indicating the paternal and maternal alleles.
    
    Args:
        snp_file (str):
            The phased SNP file.
    
    Returns:
        snp (dict):
            The phased SNP data.
    """
    snp = {}
    
    with open(snp_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            chrom = entry[0]
            loc = int(entry[2])
            p_allele, m_allele = entry[3].split('|')
            
            if chrom not in snp:
                snp[chrom] = {}
            
            snp[chrom][loc] = {}
            snp[chrom][loc]['p_allele'] = p_allele
            snp[chrom][loc]['m_allele'] = m_allele
    
    return snp


def fetch_allele_counts(pileup):
    """
    For each SNP, get the counts of each nucleotide in the overling reads.
    
    Args:
        pileup (pd.DataFrame):
            The pileup data.
    
    Returns:
        pileup_counts (pd.DataFrame):
            The nucleotide counts in each SNP pileup.
    """
    pileup_counts = []
    
    for i in range(pileup.shape[0]):
        # Get info for the SNP
        chrom = pileup['chrom'].iloc[i]
        loc = pileup['loc'].iloc[i]
        ref = pileup['ref'].iloc[i].upper()
        num = pileup['num'].iloc[i]
        
        # Initialize allele counts
        nuc = {}
        for x in ['A', 'T', 'C', 'G']:
            nuc[x] = 0
        
        for base in pileup['bases'].iloc[i]:        
            # Count number of ref alleles
            if base in ['.', ',']:
                nuc[ref] += 1
        
            # Count number of alt allele
            elif base.upper() in nuc:
                nuc[base.upper()] += 1

        pileup_counts.append([
            chrom, loc, ref, num,
            nuc['A'], nuc['T'], nuc['C'], nuc['G']])
    
    pileup_counts = pd.DataFrame(
        pileup_counts,
        columns=['chrom', 'loc', 'ref', 'num', 'A', 'T', 'C', 'G'])
    
    return pileup_counts


def compute_binomial_pvalue(pileup_counts, snp):
    """
    Compute the binomial p-value (two-tailed) for the number of 
    paternal alleles out of the number of paternal + maternal alleles.
    
    Args:
        pileup_counts (pd.DataFrame):
            The pileup allele count data.
        snp (dict):
            A dictionary of the maternal and paternal allele nucleotides
            for each phased SNP.
    
    Returns:
        res (pd.DataFrame):
            The results including p-value for each SNP.
    """
    res = []
    
    for i in range(pileup_counts.shape[0]):
        # Chrom and loc
        chrom = pileup_counts['chrom'].iloc[i]
        loc = pileup_counts['loc'].iloc[i]
        
        # Get the nucleotides of the p_allele and m_allele
        p_allele = snp[chrom][loc]['p_allele']
        m_allele = snp[chrom][loc]['m_allele']
        
        # Get allele counts
        p_count = pileup_counts[p_allele].iloc[i]
        m_count = pileup_counts[m_allele].iloc[i]
        tot_count = p_count + m_count
        
        # Get two-tailed binomial p-value
        if tot_count == 0:
            pval = 1.0
        else:
            pval = binom_test(x=p_count, n=tot_count, p=0.5)
        
        # Make allele string
        allele_string = '|'.join([p_allele, m_allele])
        
        res.append([
            chrom, loc, allele_string, p_count, m_count, pval])
    
    res = pd.DataFrame(
        res,
        columns = [
            'chrom', 'pos', 'allele', 'p_allele', 'm_allele', 'pval'])
    
    return res


def compute_phased_snp(pileup_file, snp_file, q_value_cutoff):
    """
    Compute phased SNPs.
    
    Args:
        pileup_file (str):
            The name of the pileup file from samtools (with all chroms
            combined).
        snp_file (str):
            The name of the file with phased SNPs.
        q_value_cutoff (float):
            The FDR Q-value cutoff for calling significance
    """
    # Read pileup file
    pileup = read_pileup_file(pileup_file)
    print 'read pileup'
    
    # Get allele counts from raw pileup file
    pileup_counts = fetch_allele_counts(pileup)
    print 'got allele counts'
    
    # Read SNP file
    snp = read_snp_file(snp_file)
    print 'read SNPs'
    
    # P-value
    res = compute_binomial_pvalue(pileup_counts, snp)
    print 'computed p-values'
    
    # Adjust p-value
    res['padj'] = multipletests(res['pval'], method='fdr_bh')[1]
    print 'adjusted p-values'
    
    # Remove column with raw p-value
    res.drop('pval', axis=1, inplace=True)
    print 'removed raw pvalue'
    
    # Call significant SNPs
    res['biased'] = 'No'
    res[res['padj'] < q_value_cutoff]['biased'] = 'Yes'
    print 'called significance'
    
    # Write results
    res.to_csv('phased_snp_qvalue.txt', header=True, index=False, sep='\t')


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
        '-p', '--pileup_file', required=True,
        help=('The pileup file for all chromosomes from samtools.'))
    
    required.add_argument(
        '-s', '--snp_file', required=True,
        help=('The file of phased SNPs.'))

    required.add_argument(
        '-q', '--q_value_cutoff', required=True,
        help=('The FDR Q value cutoff for significance.'))
        
    # Parse the arguments from the command-line input
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    
    # Read command-line arguments
    args = parse_command_line_args()
    
    # Compute phased SNPs
    compute_phased_snp(
        args.pileup_file,
        args.snp_file,
        float(args.q_value_cutoff))
