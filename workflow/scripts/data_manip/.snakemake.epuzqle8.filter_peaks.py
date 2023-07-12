
######## Snakemake header ########
import sys; sys.path.insert(0, "/group/flchedingrp/eth/conda/envs/sm/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X/\x00\x00\x00output/peak_data/samplePeaksDataFrame.pd.pickleq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX7\x00\x00\x00output/peak_data/samplePeaksDataFrameFiltered.pd.pickleq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12}q\x13h\x08}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17}q\x18h\x08}q\x19sbX\x07\x00\x00\x00threadsq\x1aK\x01X\t\x00\x00\x00resourcesq\x1bcsnakemake.io\nResources\nq\x1c)\x81q\x1d(K\x01K\x01e}q\x1e(h\x08}q\x1f(X\x06\x00\x00\x00_coresq K\x00N\x86q!X\x06\x00\x00\x00_nodesq"K\x01N\x86q#uh K\x01h"K\x01ubX\x03\x00\x00\x00logq$csnakemake.io\nLog\nq%)\x81q&}q\'h\x08}q(sbX\x06\x00\x00\x00configq)}q*X\x04\x00\x00\x00ruleq+X\x0c\x00\x00\x00filter_peaksq,ub.')
######## Original script #########
import pandas as pd



def subset_peaks_by_call_params(
    samp_peaks, min_length=50, min_convert=40, min_window=15
):
    """Select peaks that were called using specific call thresholds
    for length, conversion and window size.

    Args:
        samp_peaks (DataFrame): Sample dataframe with peak data
        min_length (int, optional): Min peak length. Defaults to 50.
        min_convert (int, optional): Min cytosine converion percentage. Defaults to 40.
        min_window (int, optional): Min window size. Defaults to 15.

    Returns:
        DataFrame: `samp_peaks` but only with peaks that meet thresholds
    """
    copy_df = samp_peaks.copy()

    def subset_peaks(df):
        return df[
            (df["min_len_peak"] >= min_length)
            & (df["threshold_convert"] >= min_convert)
            & (df["window"] >= min_window)
        ]

    copy_df["peaks_df"] = copy_df.apply(lambda row: subset_peaks(row.peaks_df), axis=1)

    return copy_df


def retain_unique_peaks(peaks_df):
    """Subsetting peaks by call parameters only will allow two peaks from the same read called
    with different parameters to persist. This function checks if a read has multible peak
    calls and then filters for the most stringent call. 

    Args:
        peaks_df (DataFrame): Sample dataframe with peak data

    Returns:
        DataFrame: Same dataframe but with only unique peaks
    """
    print('Processing', set(peaks_df.sample_id), 'samples')
    retained_peaks = []
    reads = set(peaks_df["read_id"])
    for each_read in reads:
        read_peaks = peaks_df.loc[peaks_df["read_id"] == each_read]
        # retain the peak that was called with the most strict parameters
        read_peaks = read_peaks.sort_values(
            by=["min_len_peak", "threshold_convert", "window"]
        )
        retained_peaks.append(read_peaks.iloc[0])
    return pd.DataFrame(retained_peaks)


if __name__ == "__main__":
    
    samp_df = pd.read_pickle(snakemake.input[0])
    subset_samp_df = subset_peaks_by_call_params(samp_df)
    subset_samp_df['peaks_df'] = subset_samp_df.apply(
        lambda row: retain_unique_peaks(row['peaks_df']),
        axis=1
    )

    subset_samp_df.to_pickle(subset_samp_df.output[0])


