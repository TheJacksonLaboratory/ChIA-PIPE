import pandas as pd


def filter_columns(d, cell_types):
    """
    """
    cell_meta = [''] + list(d.loc['Cell_type'])
    factor_meta = [''] + list(d.loc['Factor'])

    ctcf_list = [
        i for i in range(len(cell_meta)) if cell_meta[i] in cell_types and
        factor_meta[i] == 'CTCF']

    rnapII_list = [
        i for i in range(len(cell_meta)) if cell_meta[i] in cell_types and
        factor_meta[i] == 'RNAPII']
        
    d_ctcf = d[ctcf_list]
    d_rnapII = d[rnapII_list]
    
    return d_ctcf, d_rnapII 
    


if __name__ == '__main__':

    long_file = 'combined_qc_table_hg38_2017-11-20_protocol.tsv'
    situ_file = 'combined_qc_table_hg38_2018-03-15_protocol.tsv'
    cell_types = ['A549', 'GM12878', 'HUVEC', 'K562', 'NB4']
    
    long = pd.read_csv(long_file, sep='\t', index_col=0, header=None)
    situ = pd.read_csv(situ_file, sep='\t', index_col=0, header=None)
    
    long_ctcf, long_rnapII = filter_columns(long, cell_types)
    situ_ctcf, situ_rnapII = filter_columns(situ, cell_types)
    
    
    ctcf = pd.concat([long_ctcf, situ_ctcf], axis=1)
    rnapII = pd.concat([long_rnapII, situ_rnapII], axis=1)
    
    ctcf.to_csv('qc_comparison_table_by_factor_CTCF.tsv', sep='\t', header=False, index=True)
    
    
    rnapII.to_csv('qc_comparison_table_by_factor_RNAPII.tsv', sep='\t', header=False, index=True)
    