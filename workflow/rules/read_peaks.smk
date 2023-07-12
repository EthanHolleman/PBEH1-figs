

rule read_bed_files:
    conda:
        '../envs/py.yml'
    input:
        bed_dir='../resources/PBEH1-data/BED',
        sample_table='../resources/PBEH1-data/PBEH-1-samples.tsv'
    output:
        'output/peak_data/samplePeaksDataFrame.pd.pickle'
    script:'../scripts/data_manip/process_peaks.py'


rule filter_peaks:
    conda:
        '../envs/py.yml'
    input:
        'output/peak_data/samplePeaksDataFrame.pd.pickle'
    output:
        'output/peak_data/samplePeaksDataFrameFiltered.pd.pickle'
    script:'../scripts/data_manip/filter_peaks.py'


rule format_read_counts:
    conda:
        '../envs/py.yml'
    input:
        '../resources/PBEH1-data/ReadCount.tsv'
    output:
        'output/read_counts/readCountsLong.tsv'
    script:'../scripts/data_manip/parse_read_counts.py'


rule calculate_vr_rloop_frequencies:
    conda:
        '../envs/py.yml'
    input:
        read_counts='output/read_counts/readCountsLong.tsv',
        filtered_peaks= 'output/peak_data/samplePeaksDataFrameFiltered.pd.pickle'
    output:
        'test'
    script:'../scripts/data_manip/vr_loop_frequency_tables.py'





