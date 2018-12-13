###/projects/capurd/chia_pet/chia_pet_tool_2/conda/bin/python
import os
import commands
import argparse


def split_loops_into_anchors(file_name):
    """
    Split loops into anchors.
    """
    left = []
    right = []
    run = file_name.split('.')[0]
    i = -1
    
    with open(file_name, 'r') as f:
        for line in f:
            if '#' in line:
                # Skip header
                continue
            
            i += 1
            loop_id = str(i)
            entry = line.strip().split()
            
            # Left anchor
            chrom_a = entry[0]
            start_a = entry[1]
            end_a = entry[2]
            
            left.append(
                [chrom_a, start_a, end_a, loop_id, '0', '+'])
            
            # Right anchor
            chrom_b = entry[3]
            start_b = entry[4]
            end_b = entry[5]
            
            right.append(
                [chrom_b, start_b, end_b, loop_id, '0', '+'])            
    
    ## Write left anchors to BED file
    left_file = 'temp_left_anchors.bed'
    
    with open(left_file, 'w') as out:
        for row in left:
            out.write('\t'.join(row) + '\n')
    
    
    ## Write right anchors to BED file
    right_file = 'temp_right_anchors.bed'
    
    with open(right_file, 'w') as out:
        for row in right:
            out.write('\t'.join(row) + '\n')
    
    return left_file, right_file 


def bedtools_intersect(anchor_bed, peaks_file, inter_bed):
    """
    bedtools intersect.
    """
    # Remove previous overlap BED file for extra caution
    if os.path.exists(inter_bed):
        os.remove(inter_bed)
    
    command = ('module load bedtools/2.26.0;'
               ' bedtools intersect -u -a {anchor_bed} -b {peaks_file}'
               ' > {inter_bed}').format(
                    anchor_bed=anchor_bed,
                    peaks_file=peaks_file,
                    inter_bed=inter_bed)
    
    #p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    #out, err = p.communicate()
    commands.getoutput(command)
    

def annotate_loops_with_peak_support(
        loop_file, peaks_file):
    """
    Annotate loops with peak support (0, 1, 2).
    """
    left_anchor_file, right_anchor_file = \
        split_loops_into_anchors(loop_file)
    
    # Set BED file names
    left_inter_file = 'temp_left_anchors_intersect_peaks.bed'
    right_inter_file = 'temp_right_anchors_intersect_peaks.bed'
    
    # Intersect two anchors with peaks
    bedtools_intersect(left_anchor_file, peaks_file, left_inter_file)
    bedtools_intersect(right_anchor_file, peaks_file, right_inter_file)
    
    pk_supp = {}
    
    for inter_file in [left_inter_file, right_inter_file]:
        with open(inter_file, 'r') as f:
            for line in f:
                entry = line.strip().split()
                idx = int(entry[3])
                if idx not in pk_supp:
                    pk_supp[idx] = 1
                else:
                    pk_supp[idx] += 1
    
    # Read Loop file
    res = []
    with open(loop_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            res.append(entry)
    
    ## Write results
    out_file = loop_file + '.peak_annot'
    with open(out_file, 'w') as f:
        idx = -1
        for entry in res:
            idx += 1
            if idx not in pk_supp.keys():
                count = 0
            else:
                count = pk_supp[idx]
            
            entry.append(str(count))
            f.write('\t'.join(entry) + '\n')


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
		'-p', '--peak_file', required=True,
		help=('The file of peaks.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    
    annotate_loops_with_peak_support(
        args.loop_file, args.peak_file)
    