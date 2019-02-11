import argparse


def annotate_enhancer_promoter_loops(
    left_enhancers, right_enhancers, left_promoters, right_promoters):
    """
    Split loops into anchors.
    """
    loops = {}
    left_coords = {}
    right_coords = {}
    
    # Left enhancers
    with open(left_enhancers, 'r') as f:
        for line in f:
            entry = line.strip().split()
            id = entry[3]
            
            loops[id] = {}
            loops[id]['left'] = ['E']
            
            left_coords[id] = [entry[0], entry[1], entry[2]]
            
    
    # Right enhancers
    with open(right_enhancers, 'r') as f:
        for line in f:
            entry = line.strip().split()
            id = entry[3]
            
            if id in loops:
                loops[id]['right'] = ['E']
            else:
                loops[id] = {}
                loops[id]['right'] = ['E']
            
            right_coords[id] = [
                entry[0], entry[1], entry[2], entry[4], entry[5]]
            
    # Left promoters
    with open(left_promoters, 'r') as f:
        for line in f:
            entry = line.strip().split()
            id = entry[3]
            
            if id in loops:
                if 'left' in loops[id]:
                    loops[id]['left'].append('P')
                else:
                    loops[id]['left'] = ['P']
            else:
                loops[id] = {}
                loops[id]['left'] = ['P']
            
            left_coords[id] = [entry[0], entry[1], entry[2]]

    # Right promoters
    with open(right_promoters, 'r') as f:
        for line in f:
            entry = line.strip().split()
            id = entry[3]
            
            if id in loops:
                if 'right' in loops[id]:
                    loops[id]['right'].append('P')
                else:
                    loops[id]['right'] = ['P']
            else:
                loops[id] = {}
                loops[id]['right'] = ['P']
            
            right_coords[id] = [
                entry[0], entry[1], entry[2], entry[4], entry[5]]


    # Write results
    out_file = left_enhancers.replace(
        '.left_anchors.enhancers', '.enhancer_promoter')
    
    with open(out_file, 'w') as out:
        for id in sorted(loops.keys()):
            
            if 'P' in loops[id]['left'] and 'P' in loops[id]['right']:
                annot = 'PP'
            elif 'E' in loops[id]['left'] and 'P' in loops[id]['right']:
                annot = 'EP'
            elif 'P' in loops[id]['left'] and 'E' in loops[id]['right']:
                annot = 'PE'
            else:
                continue
            
            res = left_coords[id] + right_coords[id] + [annot]
            out.write('\t'.join(res) + '\n')
        

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
		'-a', '--left_enhancers', required=True,
		help=('Left anchors overlapping enhancers.'))
	
	required.add_argument(
		'-b', '--right_enhancers', required=True,
		help=('Right anchors overlapping enhancers.'))
	
	required.add_argument(
		'-c', '--left_promoters', required=True,
		help=('Left anchors overlapping promoters.'))
	
	required.add_argument(
		'-d', '--right_promoters', required=True,
		help=('Right anchors overlapping promoters.'))
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    
    annotate_enhancer_promoter_loops(
        args.left_enhancers,
        args.right_enhancers,
        args.left_promoters,
        args.right_promoters)
