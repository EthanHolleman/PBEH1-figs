# Read and process SMRF-seq called peaks (bed files). Collects all peak files
# provided as apart of PBEH-1_data directory and combines with the sample
# label table to produce a pickled pandas dataframe that is the sample label
# data combined with peak data. Additionally peaks are filtered based on their
# call parameters. By default


import pandas as pd
import re
from pathlib import Path


def stella_to_ethan_nomincature(stella_file_name):
    """Convert sample ID names assigned by Stella to their equivalent
    Ethan version (Ethan is writing this code).

    Args:
        stella_file_name (str): File name is Stella PacBio sample format

    Raises:
        Exception: Raised when no valid filename found

    Returns:
        int: Integer sample ID with no extra decoration (Ethan format)
    """
    match = re.search(r"PCB\d(\d+)", str(stella_file_name))
    if match:
        return int(match[1])
    else:
        raise Exception("No sample ID found in file string!")


def merge_bed_dir_into_sample_df(bed_dirs, samp_df):
    """Add a column into the sample dataframe that contains the path
    to directory containing all bed files for this sample.

    Args:
        bed_dirs (str): Location of all bed files for all samples
        samp_df (DataFrame): Pandas dataframe of sample information

    Returns:
        DataFrame: New sample dataframe with sample specific bed file paths
    """

    # create temp dataframe to hold sample ids and path to bed file directory
    df_list = []
    for each_path in bed_dirs:
        sample_id = stella_to_ethan_nomincature(each_path)
        df_list.append({"sample_id": sample_id, "bed_path": each_path})
    df_bed_id = pd.DataFrame(df_list)

    # merge and return dataframes
    return pd.merge(samp_df, df_bed_id, on="sample_id")


def read_PBEH_bed(filepath):
    """Parses the bed file style currently being used for PBEH1 peak calls.
    Format is bit different from standard bed to include details on parameters
    that were used for peak calling.

    Args:
        filepath (str): Path to PBEH1 bed file with peak calls

    Returns:
        DataFrame: Pandas DataFrame with labeled peak call data
    """
    bed = pd.read_csv(str(filepath), sep="\t", header=None)
    bed.columns = [
        "plasmid",
        "peak_start",
        "peak_end",
        "read_id",
        "score",
        "strand",
        "call_path",
        "UNKNOWN_VAL_0",
        "call_params",
        "call_type",
    ]
    call_params = ["threshold_convert", "window", "min_len_peak"]
    bed[call_params] = bed["call_params"].str.split(",", expand=True)

    # The call params values are left with a character in the front
    # of a numeric value. Use a regex to pull out just the numeric
    # value and replace in the dataframe
    def parse_peak_call_params(column_name):
        def get_numeric(value):
            return float(re.search(r"\d+", value)[0])

        bed[column_name] = bed.apply(lambda row: get_numeric(row[column_name]), axis=1)

        return bed

    for each_param in call_params:
        bed = parse_peak_call_params(each_param)

    # Add a length column
    bed["peak_length"] = bed["peak_end"] - bed["peak_start"]
    if "PEAK_TEMP_C" in filepath.name:
        bed["strand"] = "-"
    elif "PEAK_C" in filepath.name:
        bed["strand"] = "+"

    return bed


def collect_sample_bed_files(sample_bed_dir, sample_id):
    """Collects all bed files for a given PBEH-1 sample (a sample will have 2
    bed files for each plasmid included in that sample one for each strand) and
    combines into a single dataframe.

    Args:
        sample_bed_dir (str): Path to data directory for sample
        sample_id (int): Sample ID as an integer

    Returns:
        DataFrame: DataFrame with all peak data for a given sample ID
    """

    # there is always an extra directory called MULTIPEAKFIX
    adjusted_path = Path(sample_bed_dir).joinpath("MULTIPEAKFIX")
    frames = []
    for each_path in adjusted_path.iterdir():
        if each_path.suffix == ".bed":
            df = read_PBEH_bed(each_path)
            frames.append(df)

    sample_df = pd.concat(frames)
    sample_df["sample_id"] = sample_id

    return sample_df


def collect_peaks_by_sample(samp_bed):
    """After modifying the sample assignment dataframe with a path to each sample's
    peak data use this path to read all bed files for each sample by applying the
    `collect_sample_bed_files` function to each row of the dataframe.

    Creates a new column called `peaks_df` which is a DataFrame containing
    all peak data for a given sample.

    Args:
        samp_bed (DataFrame): Sample description dataframe with bed file paths column

    Returns:
        DataFrame: Same dataframe passed in but now with peaks_df column
    """
    samp_bed["peaks_df"] = samp_bed.apply(
        lambda row: collect_sample_bed_files(row["bed_path"], row["sample_id"]), axis=1
    )
    return samp_bed


if __name__ == "__main__":
    bed_dir = snakemake.input["bed_dir"]
    sample_table = snakemake.input["sample_table"]

    # read sample df
    samp_df = pd.read_csv(sample_table, sep="\t")
    samp_df['sample_id'] = samp_df['Sample-ID']

    bed_dirs = list(Path(bed_dir).iterdir())
    samp_bed = merge_bed_dir_into_sample_df(bed_dirs, samp_df)

    samp_peaks = collect_peaks_by_sample(samp_bed)

    samp_peaks.to_pickle(snakemake.output[0])
