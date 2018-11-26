#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import numpy as np
import argparse


def get_loop_anchor_midpoints(loop_file, min_pet_count=0):
    """
    """
    try:
        min_pet_count = int(min_pet_count)
    except:    
        min_pet_count = 0
    
    res = []
    with open(loop_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            chrom = entry[0]
            l_start = int(entry[1])
            l_end = int(entry[2])
            r_start = int(entry[4])
            r_end = int(entry[5])
            pet_count = int(entry[6])
            
            if pet_count < min_pet_count:
                continue
            
            left = int(np.mean([l_start, l_end]))
            right = int(np.mean([r_start, r_end]))
            
            res.append([chrom, left, right, pet_count])
    
    out_file = loop_file + '.anchor_mids'
    with open(out_file, 'w') as f:
        for entry in res:
            f.write('\t'.join(map(str, entry)) + '\n')



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
		help=('The file of intrachrom loops with PET count >= 2.'))
	
	# Indicate the optional arguments
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument(
		'-m', '--min_pet_count', required=False,
		help=('The minimum PET count for loops to be retained.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args
  

if __name__ == '__main__':
    
    args = parse_command_line_args()
    
    get_loop_anchor_midpoints(
        args.loop_file, args.min_pet_count)
