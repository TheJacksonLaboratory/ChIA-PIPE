#!/usr/bin/python

import os
import sys
import argparse
import pysam


def split_bam_by_species(genome_name, ctrl_genome_name):
    """
    Split the BAM file by species into one file for the primary species,
    one file for the spike-in control species, and one file for
    inter-species ligations.
    
    Args:
        genome_name (str):
            The genome name of the primary species.
        ctrl_genome_name (str):
            The genome name of the spike-in control species.
    """
    pass


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
		'-g', '--genome', required=True,
		help=('The name of the the BAM file with paired-end reads')
	required.add_argument(
		'-g', '--genome', required=True,
		help=('The name of the reference genome for the true data.')
	required.add_argument(
		'-c', '--ctrlgenome', required=True,
		help='The name of the reference genome for the spike-in control data.')
	
	# Indicate the optional arguments if any
	# None
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
	
	args = parse_command_line_args()
	split_bam_by_species(args.genome, args.ctrl_genome)
