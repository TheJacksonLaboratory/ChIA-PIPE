#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import pyBigWig
import argparse


def read_loop_file(loop_file):
    """
    """
    loops = []
    with open(loop_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            loops.append(entry)
    
    return loops


def annotate_loops_with_pileup_heights(loop_file, bigwig_file):
    """
    """
    loops = read_loop_file(loop_file)
    bw = pyBigWig.open(bigwig_file)
    
    for i in range(len(loops)):
        
        chrom_a = loops[i][0]
        start_a = int(loops[i][1])
        end_a = int(loops[i][2])
        
        chrom_b = loops[i][3]
        start_b = int(loops[i][4])
        end_b = int(loops[i][5])
        
        pileup_a = bw.stats(chrom_a, start_a, end_a, type='max')[0]
        pileup_b = bw.stats(chrom_b, start_b, end_b, type='max')[0]
        
        loops[i] += [str(int(pileup_a)), str(int(pileup_b))]
    
    
    out_file = loop_file + '.pileup_annot'
    with open(out_file, 'w') as f:
        for loop in loops:
            f.write('\t'.join(loop) + '\n')
            

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
	
	required.add_argument(
		'-b', '--bigwig_file', required=True,
		help=('The BigWig file.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args
  

if __name__ == '__main__':
    
    args = parse_command_line_args()
    
    annotate_loops_with_pileup_heights(
        args.loop_file, args.bigwig_file)
