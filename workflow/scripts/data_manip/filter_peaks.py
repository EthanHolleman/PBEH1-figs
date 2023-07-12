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


