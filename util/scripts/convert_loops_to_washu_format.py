#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import argparse
import commands


def convert_loops_to_washu_format(loop_file):
    """
    Convert loops file to format for visualization in
    WashU Epigenome Browser.
    
    args:
        loop_file (str):
            The path and name of the loop file.
    """
    ## Loop format
    # chr10	511137	513329	chr10	788654	790310	8	2
    
    ## WashU format
    # chr10	512133	512333	chr10:789382-789582,8	1	.
    # chr10	789382	789582	chr10:512133-512333,8	2	.
    
    res = []
    idx = 0
    
    # Read loop file line by line and convert to WashU format
    with open(loop_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            
            # Left anchor
            chr_a = entry[0]
            left_a = int(entry[1])
            right_a = int(entry[2])
            mid_a = int((left_a + right_a) / 2.0)
            start_a = str(mid_a - 100)
            end_a = str(mid_a + 100)
            
            coord_a = '%s:%s-%s' % (chr_a, start_a, end_a)
            
            # Right anchor
            chr_b = entry[3]
            left_b = int(entry[4])
            right_b = int(entry[5])
            mid_b = int((left_b + right_b) / 2.0)
            start_b = str(mid_b - 100)
            end_b = str(mid_b + 100)
            
            coord_b = '%s:%s-%s' % (chr_b, start_b, end_b)
            
            # PET count
            score = entry[6]
            coord_a_score = '%s,%s' % (coord_a, score)
            coord_b_score =  '%s,%s' % (coord_b, score)
            
            # Line 1 for WashU loop
            idx += 1
            line_1 = [chr_a, start_a, end_a, coord_b_score, str(idx), '.']
            res.append(line_1)
            
            # Line 2 for WashU loop
            idx += 1
            line_2 = [chr_b, start_b, end_b, coord_a_score, str(idx), '.']
            res.append(line_2)
    
    # Sort coordinates
    res = sorted(res, key=lambda x: [x[0], int(x[1])])
        
    # Write WashU file
    out_file = loop_file + '.for_WashU.mid200.txt'
    with open(out_file, 'w') as f:
        for line in res:
            f.write('\t'.join(line) + '\n')
    
    # Compress WashU file
    command_1 = ('bgzip %s' % out_file)
    commands.getoutput(command_1)
    
    # Index WashU file
    command_2 = ('tabix -p bed %s.gz' % out_file)    
    commands.getoutput(command_2)


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
    convert_loops_to_washu_format(args.loop_file)
