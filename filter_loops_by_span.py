#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import os
import argparse
import numpy as np


def filter_loops_by_span(loop_file):
    """
    Filter loops by span.
    """
    loops = []
    
    with open(loop_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            
            # Left anchor
            start_a = int(entry[1])
            end_a = int(entry[2])
            
            # Right anchor
            start_b = int(entry[4])
            end_b = int(entry[5])
            
            a = np.mean([start_a, end_a])
            b = np.mean([start_b, end_b])
            
            if abs(b - a) < 50000:
                continue
            
            loops.append(line)
    
    
    out_file = loop_file + '.min_50kb'
    
    with open(out_file, 'w') as f:
        for line in loops:
            f.write(line)

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
    
    filter_loops_by_span(args.loop_file)

