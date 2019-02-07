import glob
import numpy as np
import argparse


def split_loops_into_anchors(file_name, extend=1000):
    """
    Split loops into anchors.
    """
    left = []
    right = []
    run = file_name.split('.')[0]
    i = 0
    
    with open(file_name, 'r') as f:
        for line in f:
            if '#' in line:
                # Skip header
                continue
            
            i += 1
            loop_id = '%s_%s' % (run, i)
            entry = line.strip().split()
            
            # Scores
            pet_count = entry[6]
            n_anchor = entry[7]
            
            # Left anchor
            chrom_a = entry[0]
            start_a = int(entry[1]) - extend
            start_a = str(max(0, start_a))
            
            end_a = str(int(entry[2]) + extend)
            
            left.append(
                [chrom_a, start_a, end_a, loop_id, '0', '+',
                 pet_count, n_anchor])
            
            # Right anchor
            chrom_b = entry[3]
            start_b = int(entry[4]) - extend
            start_b = str(max(0, start_b))
            
            end_b = str(int(entry[5]) + extend)
            
            right.append(
                [chrom_b, start_b, end_b, loop_id, '0', '+',
                 pet_count, n_anchor])            
    
    ## Write left anchors to BED file
    left_file = file_name + '.L_anchors.exten_%s.bed' % extend
    
    with open(left_file, 'w') as out:
        for row in left:
            out.write('\t'.join(row) + '\n')
    
    
    ## Write right anchors to BED file
    right_file = file_name + '.R_anchors.exten_%s.bed' % extend
    
    with open(right_file, 'w') as out:
        for row in right:
            out.write('\t'.join(row) + '\n')    


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
		'-l', '--loop_file', required=True,
		help=('The file of intrachrom loops.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    
    extend = 2500
    split_loops_into_anchors(args.loop_file, extend)
