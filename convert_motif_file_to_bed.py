
def convert_motif_file_to_bed(motif_file):
    """
    """
    out_file = motif_file.replace('.txt', '.bed')
    res = []
    i = 0
    
    with open(motif_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            if entry[0] == 'ID':
                # header
                continue
            
            i += 1
            coord = entry[2]
            chrom, span = coord.split(':')
            start, end = span.split('-')
            name = 'ctcf_motif_%s' % i
            
            bed_line = [chrom, start, end, name, '0', '+']
            res.append(bed_line)
    
    with open(out_file, 'w') as f:
        for entry in res:
            f.write('\t'.join(entry) + '\n')
            

if __name__ == '__main__':
    
    motif_file = 'ctcf_motif_db_allcomp_hg18.txt'
    convert_motif_file_to_bed(motif_file)
