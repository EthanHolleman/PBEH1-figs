import pandas as pd
from Bio import SeqIO

# Prepare dataframes for generating plots about VR samples R-loop frequency, peak length etc


def subset_vr_samples(samp_peaks):
    """Select only samples with sample IDs from VR insert samples and return
    as a dataframe

    Args:
        samp_peaks (DataFrame): Dataframe with peak data

    Returns:
        DataFrame: Dataframe with only VR samples
    """

    vr_samples = set(range(1, 13))
    vr_samps = samp_peaks.loc[samp_peaks["Sample-ID"].isin(vr_samples)]
    return vr_samps


def identify_vr_location():
    anchor_region = "caaacactccctcgg"  # directly before the start of each VR
    three_prime_arm = "gaattcgtcgcagtgaccgaggcgaggagg"  # directly after each VR
    example_VR = "PBEH1-data/230513_FASTA_BED/230513_BED/FASTA/T7_init_VR_10.fa"
    vr = SeqIO.read(example_VR, "fasta")
    vr.annotations["molecule_type"] = "DNA"
    vr_start = vr.seq.find(anchor_region.upper()) + len(anchor_region)
    vr_end = vr.seq.find(three_prime_arm.upper())
    assert vr_end - vr_start == 200  # VRs are 200 bp long

    return vr_start, vr_end


def annotate_vr_intersecting(vr_samps, vr_start, vr_end):
    """Add column for each peak to indicate whether or not it intersects
    the VR insert region. Label with 1 if yes 0 if no.

    Args:
        vr_samps (DataFrame): Dataframe of VR samples with peak data

    Returns:
        DataFrame: Same dataframe but peak dataframes now annotated with
        `vr_intersecting` column
    """

    peak_dfs = []

    def is_vr_intersecting(row, vr_start, vr_end):
        if row.peak_start >= vr_start and row.peak_start < vr_end:
            return 1
        if row.peak_end >= vr_start and row.peak_end < vr_end:
            return 1
        if row.peak_start < vr_start and row.peak_end > vr_end:
            return 1

        return 0

    for i, each_row in vr_samps.iterrows():
        peak_df = each_row.peaks_df
        peak_df["vr_intersecting"] = peak_df.apply(
            lambda row: is_vr_intersecting(row, vr_start, vr_end), axis=1
        )
        peak_dfs.append(peak_df)

    vr_samps["peaks_df"] = peak_dfs

    return vr_samps


def annotate_snrpn_intersecting(vr_samps, vr_end):
    """Add a column to each peak to mark if it intersects with the first 200 bp
    of the SNRPN region downstream of the VR sequence.

    Args:
        vr_samps (DataFrame): Sample dataframe with peak data
        vr_end (int): bp location of end of VR with respect to forward primer

    Returns:
        DataFrame: Same dataframe but with additional column in peaks_df annotating SNRPN
        intersecting peaks
    """

    peaks_snrpn_dfs = []

    def is_SNRPN_intersecting(row, SNRPN_start, SNRPN_end):
        if row.peak_start >= SNRPN_start and row.peak_start < SNRPN_end:
            return 1
        
        return 0

    for i, each_row in vr_samps.iterrows():
        peak_df = each_row.peaks_df
        peak_df['SNRPN_intersecting'] = peak_df.apply(
            lambda row: is_SNRPN_intersecting(row, vr_end, vr_end+200), axis=1
        )
        peaks_snrpn_dfs.append(peak_df)

    vr_samps['peaks_df'] = peaks_snrpn_dfs

    return vr_samps


def count_peaks_by(peaks_df_vr, count_variable='vr_intersecting'):
    """Count number of peaks by a specific column. Peaks will always be grouped
    by plasmid, strand and their sample ID before counting.

    Args:
        peaks_df_vr (DataFrame): _description_
        count_variable (str, optional): _description_. Defaults to 'vr_intersecting'.

    Returns:
        _type_: _description_
    """
    count_df = peaks_df_vr.groupby(['plasmid', 'strand', 'sample_id']).count()[count_variable].reset_index()
    # add a column that has the plasmid ID as upper case, use this for merge bc
    # this is how plasmid names are formatted in read count dataframe
    count_df['plasmid_upper'] = count_df.apply(
        lambda row: row.plasmid.upper(),
        axis=1
    )
    return count_df

def subset_vr_intersecting(peaks_df):
    return peaks_df.loc[peaks_df['vr_intersecting'] == 1]


def subset_snrpn_intersecting(peaks_df):
    return peaks_df.loc[peaks_df['SNRPN_intersecting'] == 1]


def merge_total_reads(vr_peak_count, mod_read_df):
    """Merge total read statistics into peak count dataframes in order to be able to
    calculate R-loop frequencies.

    Args:
        vr_peak_count (DataFrame): Peak count dataframe
        mod_read_df (DataFrame): Read count dataframe in long format. See `parse_read_counts.py` script

    Returns:
        DataFrame: Peak count dataframe merged with read count
    """
    
    return vr_peak_count.merge(
        mod_read_df, left_on=['plasmid_upper', 'strand', 'sample_id'], 
        right_on=['Gene', 'Strand', 'BC']
    )


def calculate_vr_rloop_freq(peak_count_df, count_column, freq_column_name):
    
    peak_count_df[freq_column_name] = peak_count_df.apply(
        lambda row: round(row[count_column] / row.Reads, 3),
        axis=1
    )
    return peak_count_df





if __name__ == '__main__':

    samp_peaks = pd.read_pickle(snakemake.input['samp_peaks'])
    vr_samps = subset_vr_samples(samp_peaks)

    # identify VR start and end positions
    vr_start, vr_end = 

    # annotate peaks that intersect VRs
    vr_samps = annotate_vr_intersecting(vr_samps, vr_start, vr_end)
    vr_samps = annotate_snrpn_intersecting(vr_samps, vr_end)

    # create additional columns that contain only peaks that are either
    # VR or start of SNRPN intersecting
    vr_samps['peaks_df_vr'] = vr_samps.apply(lambda row: subset_vr_intersecting(row.peaks_df), axis=1)
    vr_samps['SNRPN_peaks'] = vr_samps.apply(lambda row: subset_snrpn_intersecting(row.SNRPN_peaks), axis=1)

    # Create new column that contains peak counts for either VR or SNRPN intersecting peaks
    vr_samps['vr_peak_count'] = vr_samps.apply(lambda row: count_peaks_by(row.peaks_df_vr), axis=1)
    vr_samps['snrpn_peak_count'] = vr_samps.apply(lambda row: count_peaks_by(row.SNRPN_peaks, 'SNRPN_intersecting'), axis=1)

    # merge in total read cound data into peak counts so can next caculate
    # peak frequency
    vr_samps['vr_peak_count'] = vr_samps.apply(
        lambda row: merge_total_reads(row.vr_peak_count, mod_read_df),
        axis=1
    )

    vr_samps['snrpn_peak_count'] = vr_samps.apply(
        lambda row: merge_total_reads(row.snrpn_peak_count, mod_read_df),
        axis=1
    )

    # Calculate R-loop frequency for out regions of interest (VR and SNRPN)
    vr_samps['vr_peak_count'] = vr_samps.apply(
        lambda row: calculate_vr_rloop_freq(row.vr_peak_count, 'vr_intersecting', 'vr_rloop_freq'),
        axis=1
    )

    vr_samps['snrpn_peak_count'] = vr_samps.apply(
        lambda row: calculate_vr_rloop_freq(row.snrpn_peak_count, 'SNRPN_intersecting', 'SNRPN_rloop_freq'),
        axis=1
    )







