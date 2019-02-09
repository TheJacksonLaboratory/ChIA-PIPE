import os
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.SeqRecord import SeqRecord
import regex
import argparse
import gzip


def parse_read_pair_for_linker(r1, r2, linker, n_mismatch, min_tag_len):
    """
    Parse a read pair (R1, R2) for the linker sequence.
    
    Args:
        r1 (Bio.SeqRecord.SeqRecord):
            R1 read as a Biopython record.
        r2 (Bio.SeqRecord.SeqRecord):
            R2 read as a Biopython record.
        linker (str):
            The linker sequence.
        n_mismatch (int):
            The number of mismatches allowed when searching for
            the linker.
        min_tag_len (int):
            The minimum tag length when deciding if a read pair is
            single tag or PET. 
            
    Returns:
        r1 (Bio.SeqRecord.SeqRecord):
            The parsed R1 read, where the sequence is now just the
            genomic tag.
        r2 (Bio.SeqRecord.SeqRecord):
            The parsed R2 read, where the sequence is now just the
            genomic tag.
        n_link (int):
            The number of reads in which a linker was found: 0, 1, or 2.
        n_tag (int):
            The number of usable genomic tags (length >=18 bp): 0, 1, or 2.
    """
    # Set minimum tag length
    min_tag_len = 18
    
    # Initialize number of linkers to zero
    n_link = 0
    
    # Initialize number of tags to zero
    n_tag = 0
    
    # Check R1 for linker with exact match
    exact_r1 = regex.findall(linker, str(r1.seq), overlapped=True)
    
    # Check R1 for linker with 1 mismatch
    mm_r1 = regex.findall(
        "(%s){s<=%s}" % (linker, n_mismatch),
        str(r1.seq),
        overlapped=True)
    
    # Check R2 for linker with exact match
    exact_r2 = regex.findall(linker, str(r2.seq), overlapped=True)
    
    # Check R2 for linker with 1 mismatch
    mm_r2 = regex.findall(
        "(%s){s<=%s}" % (linker, n_mismatch),
        str(r2.seq),
        overlapped=True)
    
    # Initialize flags
    r1_split = None    
    r2_split = None
    
    # Case 1A: R1 had an exact-match linker
    if len(exact_r1) > 0:
        n_link += 1
        r1_split = r1.seq.split(linker)
    
    # Case 1B: R1 had a 1bp-mismatch linker
    elif len(mm_r1) > 0:
        n_link += 1
        r1_split = r1.seq.split(mm_r1[0])
    
    # R1 can be split (Case 1A or 1B True)
    if r1_split:
        
        # Tag is the 5' subsequence
        new_r1_seq = r1_split[0]
        
        # Post-process to avoid error
        # Sequence must have len > 0 for Biopython
        if len(new_r1_seq) == 0:
            new_r1_seq = Seq('A')
          
        new_r1_qual = \
            r1.letter_annotations['phred_quality'][:len(new_r1_seq)]
        
        r1.letter_annotations = {}
        r1.seq = new_r1_seq
        r1.letter_annotations['phred_quality'] = new_r1_qual
    
    
    # Case 2A: R2 had an exact-match linker
    if len(exact_r2) > 0:
        n_link += 1
        r2_split = r2.seq.split(linker)
    
    # Case 2B: R1 had a 1bp-mismatch linker
    elif len(mm_r2) > 0:
        n_link += 1
        r2_split = r2.seq.split(mm_r2[0])
    
    # R2 can be split (Case 2A or 2B True)
    if r2_split:
        
        # Tag is the 5' subsequence
        new_r2_seq = r2_split[0]
        
        # Post-process to avoid error
        # Sequence must have len > 0 for Biopython
        if len(new_r2_seq) == 0:
            new_r2_seq = Seq('A')
        
        new_r2_qual = \
            r2.letter_annotations['phred_quality'][-len(new_r2_seq):]
        
        r2.letter_annotations = {}
        r2.seq = new_r2_seq
        r2.letter_annotations['phred_quality'] = new_r2_qual
    
    
    # Test length of genomic tag 1
    if len(r1.seq) >= min_tag_len:
        n_tag += 1
    # Tes length of genomic tag 2
    if len(r2.seq) >= min_tag_len:
        n_tag += 1
                
    return r1, r2, n_link, n_tag


def standardize_read_name(r):
    """
    Standardize the read name.
    
    Args:
        r (Bio.SeqRecord.SeqRecord):
            The read as a Biopython record.
            
    Returns:
        r (Bio.SeqRecord.SeqRecord):
            The read with a standardized name. 
    """
    old_name_parts = r.name.split('.')
    
    read_num = old_name_parts[1]
    mate_num = old_name_parts[2]
    
    if mate_num == '1':
        mate_sect = '90:130'
    else:
        mate_sect = '94:129' 
      
    generic_sect = '897:1101'
    
    read_name = ':'.join(
        [read_num, mate_sect, generic_sect, read_num, read_num])
    
    r.description = read_name
    r.name = read_name
    r.id = read_name
                
    return r


def open_output_files(run):
    """
    Open output files in appending mode.
    """
    none_file = '%s.none.fastq' % run
    one_tag_file = '%s.singlelinker.single.fastq' % run
    two_tag_file = '%s.singlelinker.paired.fastq' % run
    conflict_file = '%s.conflict.fastq' % run
    
    if os.path.exists(none_file):
        os.remove(none_file)
    
    if os.path.exists(one_tag_file):
        os.remove(one_tag_file)
        
    if os.path.exists(two_tag_file):
        os.remove(two_tag_file)
        
    if os.path.exists(conflict_file):
        os.remove(conflict_file)
    
    a_none = open(none_file, 'a') 
    a_one_tag = open(one_tag_file, 'a') 
    a_two_tag = open(two_tag_file, 'a') 
    a_conflict = open(conflict_file, 'a')
    
    return a_none, a_one_tag, a_two_tag, a_conflict


def close_output_files(a_none, a_one_tag, a_two_tag, a_conflict):
    """
    Close output files.
    """
    a_none.close()
    a_one_tag.close()
    a_two_tag.close()
    a_conflict.close()


def filter_hichip_linker(
    r1_file, r2_file, linker, n_mismatch, min_tag_len, run):
    """
    Filter read pairs for HiChIP linker and partition into categories.
    
    Args:
        r1_file (str):
            The path/name of the R1 FASTQ file.
        r2_file (st)::
            The path/name of the R2 FASTQ file.
        linker (str):
            The HiChIP linker sequence (e.g., GATCGATC)
        n_mismatch (int):
            The number of mismatches allowed when searching for
            the linker.
        min_tag_len (int):
            The minimum tag length when deciding if a read pair is
            single tag or PET. 
        run (str):
            The sequencing run ID.
        out_dir (str):
            The output directory to which to write the results.
    """   
    # Get R1 reads
    #r1_list = list(SeqIO.parse(r1_file, "fastq"))
    r1_iter = FastqPhredIterator(gzip.open(r1_file, 'rt'))
    
    # Get R2 reads
    #r2_list = list(SeqIO.parse(r2_file, "fastq"))
    r2_iter = FastqPhredIterator(gzip.open(r2_file, 'rt'))
    
    # Initialize results dictionary for partitioning reads
    
    # Open files
    a_none, a_one_tag, a_two_tag, a_conflict = open_output_files(run)
    
    # Iterate over read pairs
    for r1, r2 in itertools.izip(r1_iter, r2_iter):
        
        ## Parse the read pair for the linker
        # The sequence of each read is now just the genomic tag
        # The number of linkers and number of tags is tracked
        r1, r2, n_link, n_tag = parse_read_pair_for_linker(
            r1, r2, linker, n_mismatch, min_tag_len)
        r1 = standardize_read_name(r1)
        r2 = standardize_read_name(r2)
        res = [r1, r2]
        
        ## Categorize read pair appropriately
        if n_link == 0:
            # No linker
            SeqIO.write(res, a_none, 'fastq')
        elif n_link == 2:
            # Two or more linkers -> conflict
            SeqIO.write(res, a_conflict, 'fastq')
        elif n_link == 1:
            # Single linker
            if n_tag == 0:
                # No tags -> conflict
                SeqIO.write(res, a_conflict, 'fastq')
            elif n_tag == 1:
                # One tag
                SeqIO.write(res, a_one_tag, 'fastq')
            elif n_tag == 2:
                # Two tags
                SeqIO.write(res, a_two_tag, 'fastq')
    
    close_output_files(a_none, a_one_tag, a_two_tag, a_conflict)


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
		'-a', '--r1_file', required=True,
		help=('The path/filename of the R1 FASTQ file.'))
	
	required.add_argument(
		'-b', '--r2_file', required=True,
		help=('The path/filename of the R2 FASTQ file.'))
	
	required.add_argument(
		'-r', '--run', required=True,
		help=('The ID of the run.'))
    
	required.add_argument(
		'-l', '--linker', required=True,
		help=('The linker sequence (e.g., GATCGATC).'))
	
	required.add_argument(
		'-m', '--min_tag_len', required=True,
		help=('The minimum tag length for mapping (e.g., 18bp).'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    n_mismatch = 2
    
    # Perform the HiChIP linker filtering
    filter_hichip_linker(
        args.r1_file,
        args.r2_file, 
        args.linker,
        n_mismatch,
        args.min_tag_len,
        args.run)
    