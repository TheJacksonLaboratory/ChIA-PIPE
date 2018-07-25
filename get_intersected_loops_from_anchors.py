

def parse_left_anchor_intersects(left_file):
    """
    Parse file of left-anchor intersects.
    """
    loop_coords = {}
    intersects = {}
    
    with open(left_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            
            # Add Rep1 Loop ID to loop_coords dict
            loop_id_1 = entry[3]
            if loop_id_1 not in loop_coords: 
                loop_coords[loop_id_1] = {}
                loop_coords[loop_id_1]['chrom'] = entry[0] 
                loop_coords[loop_id_1]['start'] = entry[1]
                loop_coords[loop_id_1]['end'] = entry[2]
                loop_coords[loop_id_1]['q'] = entry[6]
                loop_coords[loop_id_1]['pet_count'] = entry[7]
                loop_coords[loop_id_1]['other'] = entry[8]
            
            # Add Rep2 Loop ID to loop_coords dict
            loop_id_2 = entry[12]
            if loop_id_2 not in loop_coords: 
                loop_coords[loop_id_2] = {}
                loop_coords[loop_id_2]['chrom'] = entry[9]
                loop_coords[loop_id_2]['start'] = entry[10]
                loop_coords[loop_id_2]['end'] = entry[11]
                loop_coords[loop_id_2]['q'] = entry[15]
                loop_coords[loop_id_2]['pet_count'] = entry[16]
                loop_coords[loop_id_2]['other'] = entry[17]
            
            # Add the pair of loop IDs to the intersects dict
            if loop_id_1 not in intersects:
                intersects[loop_id_1] = {}
            
            intersects[loop_id_1][loop_id_2] = [1, 0]
    
    return loop_coords, intersects


def parse_right_anchor_intersects(right_file, loop_coords, intersects):
    """
    Parse file of right-anchor intersects.
    """
    with open(right_file, 'r') as f:
        for line in f:
            entry = line.strip().split()
            
            # Fix end position for Rep1 Loop ID
            loop_id_1 = entry[3]
            if loop_id_1 in loop_coords: 
                loop_coords[loop_id_1]['end'] = entry[2]
            
            # Fix end position for Rep2 Loop ID
            loop_id_2 = entry[12]
            if loop_id_2 in loop_coords: 
                loop_coords[loop_id_2]['end'] = entry[11]
            
            # Record if right anchor intersects as well
            if loop_id_1 in intersects:
                if loop_id_2 in intersects[loop_id_1]:
                    intersects[loop_id_1][loop_id_2] = [1, 1]
    
    return loop_coords, intersects


def create_intersected_loops_row(loop_coords, loop_id_1, loop_id_2):
    """
    Create intersected loops row.
    """
    row = [
        loop_coords[loop_id_1]['chrom'],
        loop_coords[loop_id_1]['start'],
        loop_coords[loop_id_1]['end'],
        loop_id_1,
        '0',
        '+',
        loop_coords[loop_id_1]['q'],
        loop_coords[loop_id_1]['pet_count'],
        loop_coords[loop_id_1]['other'],
        loop_coords[loop_id_2]['chrom'],
        loop_coords[loop_id_2]['start'],
        loop_coords[loop_id_2]['end'],
        loop_id_2,
        '0',
        '+',
        loop_coords[loop_id_2]['q'],
        loop_coords[loop_id_2]['pet_count'],
        loop_coords[loop_id_2]['other']]
    
    return row


def write_intersected_loops(loop_coords, intersects, out_file):
    """
    Write intersected loops.
    """
    res = []
    rep1_ids = []
    rep2_ids = []
    
    for loop_id_1 in intersects.keys():
        for loop_id_2 in intersects[loop_id_1].keys():
            
            if intersects[loop_id_1][loop_id_2] == [1, 1]:
                
                row = create_intersected_loops_row(
                    loop_coords, loop_id_1, loop_id_2)
                    
                res.append(row)
                rep1_ids.append(loop_id_1)
                rep2_ids.append(loop_id_2)
    
    # Write loop intersects "wa wb"
    with open(out_file, 'w') as f:
        for row in res:
            f.write('\t'.join(row) + '\n')
    
    # Determine minimum intersect size "u"
    rep1_ids = set(rep1_ids)
    rep2_ids = set(rep2_ids)
    u = min(len(rep1_ids), len(rep2_ids))
    print u

    
def get_intersected_loops_from_anchors(left_file, right_file):
    """
    Get intersected loops from anchors.
    """
    loop_coords, intersects = \
        parse_left_anchor_intersects(left_file)
    
    loop_coords, intersects = \
        parse_right_anchor_intersects(
            right_file, loop_coords, intersects)
    
    
    out_file = left_file.replace('left_anchors', 'loops')
    #out_file = 'loops_head'
    write_intersected_loops(loop_coords, intersects, out_file)
    

if __name__ == '__main__':
    
    left_file = 'intersect_left_anchors_LHG0010H_LHG0014H_wa_wb.txt'
    right_file = 'intersect_right_anchors_LHG0010H_LHG0014H_wa_wb.txt'
    
    get_intersected_loops_from_anchors(left_file, right_file)
    