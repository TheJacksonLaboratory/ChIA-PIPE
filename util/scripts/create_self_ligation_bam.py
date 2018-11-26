#!/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import argparse
import pysam


def create_self_ligation_bam(in_bam):
    """
    Create self-ligation BAM file.
    """
    # Initialize input BAM
    b = pysam.AlignmentFile(in_bam, 'rb')
    
    # Initialize output BAM (self-ligation reads)
    out_bam = in_bam.replace('.bam', '.self-ligation.bam')
    o = pysam.AlignmentFile(out_bam, 'wb', template=b)
    
    # Iterate over reads in input BAM
    for read in b.fetch():
        if read.rname != read.rnext:
            # Interchromosomal pair
            continue
        
        if abs(int(read.pos) - int(read.pnext)) > 8000:
            # Span > 8 kb
            continue
        
        # Write self-ligation
        #print 'Y'
        o.write(read)

    # Close files
    b.close()
    o.close()


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
		'-b', '--bam_file', required=True,
		help=('The input BAM file.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    create_self_ligation_bam(args.bam_file)
    