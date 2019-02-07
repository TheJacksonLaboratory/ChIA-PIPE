import os
import sys
import argparse
import pysam


def adjust_loop_anchor_using_bam(IN_BAM, IN_LOOP):
    """
    Adjust loop-anchor position using pileup max from BAM file.
    """
    # Open BAM file
    samfile = pysam.Samfile(IN_BAM, "rb")
    
    # Open loop file (original)
    FLOOP = open(IN_LOOP, "r")
    
    # Open output loop file (adjusted)
    FOUT = open(IN_LOOP + ".adjusted_from_bam", "w")
    
    while True:
        line = FLOOP.readline()
        if not line:
            break
        word = line.rstrip("\n").split("\t")
        
        # Parse line in loop file
        CZ1 = word[0]
        SZ1 = int(word[1])
        EZ1 = int(word[2])
        CZ2 = word[3]
        SZ2 = int(word[4])
        EZ2 = int(word[5])
        PET = word[6]	
        
        ## Left anchor
        
        # Initialize variables for getting max pileup from BAM file
        AZ1S = 0
        AZ1E = 0
        AZ1V = 0

        DZ1S = {}
        DZ1E = {}
        DZ1V = {}
        DZ1LEN = {}
        DZ1MID = {}
        MIDZ1 = 0
        
        # Find position of max pileup in left anchor
        for pileupcolumn in samfile.pileup(CZ1, SZ1, EZ1):
            
            MIDZ1=int((SZ1 + EZ1) / 2)
            
            if AZ1S == 0 and AZ1E == 0 and AZ1V == 0:
                AZ1S = pileupcolumn.pos
                AZ1E = pileupcolumn.pos
                AZ1V = pileupcolumn.n
                DZ1S[AZ1S] = AZ1S
                
                DZ1E[AZ1S] = AZ1E
                DZ1V[AZ1S] = AZ1V
                DZ1LEN[AZ1S] = AZ1E - AZ1S
                DZ1MID[AZ1S] = int((AZ1E + AZ1S) / 2)	
            
            elif pileupcolumn.pos == (AZ1E + 1) and pileupcolumn.n == AZ1V:
                AZ1E = pileupcolumn.pos
            
            else:
                DZ1S[AZ1S] = AZ1S
                DZ1E[AZ1S] = AZ1E
                DZ1V[AZ1S] = AZ1V
                DZ1LEN[AZ1S] = AZ1E - AZ1S
                DZ1MID[AZ1S] = int((AZ1E + AZ1S) / 2)
                
                AZ1S = pileupcolumn.pos
                AZ1E = pileupcolumn.pos
                AZ1V = pileupcolumn.n
        
        if DZ1V == {}:
            KMZ1 = MIDZ1
            DZ1V[SZ1] = 0
            DZ1S[SZ1] = SZ1
            DZ1E[SZ1] = EZ1
            mZ1a = SZ1
        
        else:
            KMZ1 = max(DZ1V, key=DZ1V.get)
            RZ1DIS = {}
            
            for K in DZ1V:
                if DZ1V[K] == DZ1V[KMZ1]:
                    RZ1DIS[K] = abs(MIDZ1 - DZ1MID[K])
            
            mZ1a = min(RZ1DIS, key=RZ1DIS.get)
    
        if DZ1E[mZ1a] == DZ1S[mZ1a]:
            DZ1E[mZ1a] = DZ1E[mZ1a] + 1
        
        
        ## Right anchor
        
        # Initialize variables for getting max pileup from BAM file
        AZ2S = 0
        AZ2E = 0
        AZ2V = 0

        DZ2S = {}
        DZ2E = {}
        DZ2V = {}
        
        DZ2LEN = {}
        DZ2MID = {}
        MIDZ2 = 0
        
        # Find position of max pileup in right anchor
        for pileupcolumn in samfile.pileup(CZ2, SZ2, EZ2):
            
            MIDZ2 = int((SZ2 + EZ2) / 2)
            
            if AZ2S == 0 and AZ2E == 0 and AZ2V == 0:
                AZ2S = pileupcolumn.pos
                AZ2E = pileupcolumn.pos
                AZ2V = pileupcolumn.n
                DZ2S[AZ2S] = AZ2S
                            
                DZ2E[AZ2S] = AZ2E
                DZ2V[AZ2S] = AZ2V
                DZ2LEN[AZ2S] = AZ2E - AZ2S
                DZ2MID[AZ2S] = int((AZ2E + AZ2S) / 2)	
            
            elif pileupcolumn.pos == (AZ2E + 1) and pileupcolumn.n == AZ2V:
                AZ2E = pileupcolumn.pos
            
            else:
                DZ2S[AZ2S] = AZ2S
                DZ2E[AZ2S] = AZ2E
                DZ2V[AZ2S] = AZ2V
                DZ2LEN[AZ2S] = AZ2E - AZ2S
                DZ2MID[AZ2S] = int((AZ2E + AZ2S) / 2)
                
                AZ2S = pileupcolumn.pos
                AZ2E = pileupcolumn.pos
                AZ2V = pileupcolumn.n
        
        if DZ2V == {}:
            KMZ2 = MIDZ2
            DZ2V[SZ2] = 0
            DZ2S[SZ2] = SZ2
            DZ2E[SZ2] = EZ2
            mZ2a = SZ2
        
        else:
            KMZ2 = max(DZ2V, key=DZ2V.get)
            RZ2DIS = {}
            
            for K in DZ2V:
                if DZ2V[K] == DZ2V[KMZ2]:
                    RZ2DIS[K] = abs(MIDZ2 - DZ2MID[K])
            
            mZ2a = min(RZ2DIS, key=RZ2DIS.get)
    
        if DZ2E[mZ2a] == DZ2S[mZ2a]:
            DZ2E[mZ2a] = DZ2E[mZ2a] + 1
        
        
        ## Write adjusted loop-anchor positions
        new_anchor_line = '\t'.join([
            CZ1, str(DZ1S[mZ1a]), str(DZ1E[mZ1a]),
            CZ2, str(DZ2S[mZ2a]), str(DZ2E[mZ2a]),
            str(PET)
        ])
            
        FOUT.write(new_anchor_line + '\n')
    
    # Close files
    samfile.close()
    FLOOP.close()
    FOUT.close()


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
	
	required.add_argument(
		'-l', '--loop_file', required=True,
		help=('The input loop file.'))
	
	# Parse the arguments from the command-line input
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
    
    args = parse_command_line_args()
    adjust_loop_anchor_using_bam(args.bam_file, args.loop_file)
