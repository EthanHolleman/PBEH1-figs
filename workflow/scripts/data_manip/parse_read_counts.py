import pandas as pd



def to_long_read_count(read_count_df):
    """Convert read_count table into long format so later it will be easier to merge
    with other similary formated dataframes.

    Args:
        read_count_df (DataFrame): Dataframe containing read count information for all samples
        and all plasmids.

    Returns:
        DataFrame: Long formatted version of the same dataframe with only cirtical
        columns retained.
    """
    read_counts = []

    for i, row in read_count_df.iterrows():
        row_dict_plus_strand = {
            'Gene': row.Gene,
            'BC': row.BC,
            'Strand': '+',
            'Reads': row.READ_C
        }
        row_dict_minus_strand = {
            'Gene': row.Gene,
            'BC': row.BC,
            'Strand': '-',
            'Reads': row.READ_TEMP_C
        }
        read_counts.append(row_dict_plus_strand)
        read_counts.append(row_dict_minus_strand)

    mod_read_df = pd.DataFrame(read_counts) 
    
    return mod_read_df



if __name__ == '__main__':

    read_count = pd.read_csv(snakemake.input['read_count'], sep='\t')
    read_count_long = to_long_read_count(read_count)
    read_count_long.to_csv(snakemake.output[0], sep='\t', index=False)
