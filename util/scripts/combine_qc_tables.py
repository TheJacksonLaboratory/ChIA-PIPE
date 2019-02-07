import os
import argparse
import pandas as pd
import datetime
import random
from string import ascii_lowercase


def combine_qc_matrices(library_list):
    """
    Combines the QC matrices from multiple ChIA-PET libraries.
    
    Args:
        library_list (str):
            The name of the file containing the list of libraries
            for which to combine the QC matrices.
    """
    libs = []
    with open(library_list, 'r') as f:
        for line in f:
            entry = line.strip().split()[0]
            libs.append(entry)
    
    # Initialize the data frame with the QC matrix
    # of the first library
    matrix_file = '%s/%s.final_stats.tsv' % (libs[0], libs[0])
    d = pd.read_csv(
        matrix_file, header=None, sep='\t', index_col=0,names=[libs[0]])
    
    # Load the QC matrices for the subsequent libraries one-by-one
    for i in range(1, len(libs)):
        matrix_file = '%s/%s.final_stats.tsv' % (libs[i], libs[i])
        mat = pd.read_csv(
            matrix_file, header=None, sep='\t', index_col=0,names=[libs[i]])

        # Concatenate to the data frame the QC matrix
        # for each additional library
        d = pd.concat([d, mat], axis=1, join_axes=[d.index])
    
    # Get date and letters tag
    date = str(datetime.date.today())
    #letters = ''
    #idxs = random.sample(range(25), 2)
    #for idx in idxs:
    #    letters += ascii_lowercase[idx]
    
    # Get reference genome
    genome = d.loc['Reference_genome'][d.columns[0]]
    if '/' in genome:
        genome = genome.split('/')[-1]
        if '.' in genome:
            genome = genome.split('.')[0]
    
    # Write
    out_file = "combined_qc_table_%s_%s.tsv" % (genome, date)
    d.to_csv(out_file, sep="\t", header=False, index=True)


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
		'-l', '--library_list', required=True,
		help=('A file of the library IDs for which to combine the QC '
		      'matrices (one ID per line).'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    # Parse command-line arguments
    args = parse_command_line_args()
    
    # Perform the read-level haplotyping
    combine_qc_matrices(args.library_list)
    