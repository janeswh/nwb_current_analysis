import pandas as pd
import matplotlib.pyplot as plt
import pynwb
import os
import numpy as np
from neo.io import IgorIO
from pynwb import NWBHDF5IO
import elephant

# gets sweep info for all cells
sweep_info = pd.read_csv('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables/second_dataset_sweep_info.csv', index_col=0)

# reads in NWB file
io = NWBHDF5IO('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data/JH20210922_c3.nwb', 'r', load_namespaces=True)
nwbfile = io.read()

# # gets the IC timeseries data from the first sweep
# nwbfile.sweep_table.series[0].data[:] 

# turn sweep_table into a pandas df
nwb_df = pd.DataFrame(nwbfile.sweep_table[:]) 

# # defines list of sweeps to pull from sweep table, then extracts into subset
# select_sweeps = list(range(11, 21)) 
# subset_nwb_df = nwb_df[nwb_df['sweep_number'].isin(select_sweeps)] 

# finds the series indices for VoltageClampSeries data
b = [series[0] for series in nwb_df['series']]

vc_indices = []
for index, item in enumerate(b):
    print(type(b[index]) == pynwb.icephys.VoltageClampSeries, index)
    if type(b[index]) == pynwb.icephys.VoltageClampSeries:
        vc_indices.append(index)

# puts the VC series data into another df
vc_sweeps = nwb_df.iloc[vc_indices]

# extracts VC series data into df
raw_series = []

for i in vc_indices:
    raw_series.append(b[i].data[:])

# columns are sweeps, rows are time (*fs)
raw_df = pd.DataFrame(raw_series).T
raw_df.columns = vc_sweeps.sweep_number




# gets sweep info for one cell, drops empty values
file = 'JH20210922_c3'
file_split = file.split('.')
cell_name = file_split[0]

cell_sweep_info = pd.DataFrame(sweep_info.loc[cell_name]).dropna()

# pulls the non-escaped sweep numbers out
# if (cell_sweep_info.loc['Escaped sweeps'] == 'Yes').any:
if 'Escaped sweeps' in cell_sweep_info.index:
    nonesc_range = cell_sweep_info.loc['Non-esc sweeps'][0].split(',')
    nonesc_sweeps = []
    for i in range(len(nonesc_range)):
        if '-' in nonesc_range[i]:
            r_start = int(nonesc_range[i].split('-')[0])
            r_end = int(nonesc_range[i].split('-')[1])+1
            all_sweeps = list(range(r_start, r_end))
            nonesc_sweeps.extend(all_sweeps)
        else:            
            nonesc_sweeps.append(int(nonesc_range[i]))

# pulls out depol sweep numbers
if 'Depol sweeps' in cell_sweep_info.index:
    depol_range = cell_sweep_info.loc['Depol sweeps'][0].split(',')
    depol_sweeps = []
    for i in range(len(depol_range)):
        if '-' in depol_range[i]:
            r_start = int(depol_range[i].split('-')[0])
            r_end = int(depol_range[i].split('-')[1])+1
            all_sweeps = list(range(r_start, r_end))
            depol_sweeps.extend(all_sweeps)
        else:            
            depol_sweeps.append(int(nonesc_range[i]))



# if applicable, drops depol sweeps and  esc sweeps
if 'depol_sweeps' in globals():
    raw_df = raw_df.drop(columns=depol_sweeps, axis=1)

if 'nonesc_sweeps' in globals():
    raw_df = raw_df.filter(nonesc_sweeps)



stim_sweep_info = cell_sweep_info.filter(like="%", axis=0)

# create dict with stim name as keys, VC data as values
sweeps_dict = {}

# define sweep ranges for each stim set present
for i in range(len(stim_sweep_info.index)):
    stim_name = stim_sweep_info.index[i]
    stim_range = stim_sweep_info.iloc[i].str.split(',')[0]
    stim_sweeps = []
    # if len(stim_range) == 1:
    #     if '-' not in stim_range[0]:
    #         stim_sweeps.append(int(stim_range[0]))
    #     else:
    #         r_start = int(stim_range[0].split('-')[0])
    #         r_end = int(stim_range[0].split('-')[1])+1
    #         all_sweeps = list(range(r_start, r_end))
    #         stim_sweeps.extend(all_sweeps)
    # else:
    for j in range(len(stim_range)):
        if '-' in stim_range[j]:
            r_start = int(stim_range[j].split('-')[0])
            r_end = int(stim_range[j].split('-')[1])+1
            all_sweeps = list(range(r_start, r_end))
            stim_sweeps.extend(all_sweeps)
        else:            
            stim_sweeps.append(int(stim_range[j][0]))

    stim_sweeps_VC = raw_df[raw_df.columns.intersection(set(stim_sweeps))]
    sweeps_dict[stim_name] = stim_sweeps_VC

# drop keys with empty dataframes
sweeps_dict = {k:v for (k,v) in sweeps_dict.items() if not v.empty}

# check that dict code works for different cells with different depol/esc ranges/types




''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500  # Hz
stim_time = 520     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
tp_start = 5        # ms, time of start of test pulse
vm_jump = 10        # mV, test pulse voltage jump
pre_tp = 3          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1      # scaler for making plots in pA
fs = 25             # kHz, the sampling frequency


stim_conditions = list(sweeps_dict)
genotype = cell_sweep_info.loc['Genotype']

''' Functions '''

def filter_traces(traces, fs, lowpass_freq):
    '''
    add filtered traces attrbute to data object
    lowpass_freq: frequency in Hz to pass to elephant filter
    '''
    traces_filtered = elephant.signal_processing.butter(traces.T, lowpass_freq=lowpass_freq, fs=fs * 1000)
    traces_filtered = pd.DataFrame(traces_filtered).T

    return traces_filtered


def calculate_mean_baseline(data, fs, baseline_start, baseline_end=6000):
    '''
    Find the mean baseline in a given time series, defined as the last 3s of the sweep
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    baseline_start: int or float
        The time in ms when baseline starts.
    baseline_end: int or float
        The length of the sweep and the time when baseline ends.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = baseline_start * fs
    stop = baseline_end * fs
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def calculate_std_baseline(data, fs, baseline_start, baseline_end=6000):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    baseline_start: int or float
        The time in ms when baseline starts.
    baseline_end: int or float
        The length of the sweep and the time when baseline ends.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = baseline_start * fs
    stop = baseline_end * fs
    window = data.iloc[start:stop]
    std = window.std()

    return std


def calculate_epsc_peak(data, baseline, fs, stim_time, post_stim=100, polarity='-', index=False):
    '''
    Find the peak EPSC value for a pandas.Series or for each sweep (column) of
    a pandas.DataFrame. This finds the absolute peak value of mean baseline
    subtracted data.

    Parameters:
    -----------
    
    data: pandas.Series or pandas.DataFrame
        Time series data with stimulated synaptic response triggered at the
        same time for each sweep.
    baseline: scalar or pandas.Series
        Mean baseline values used to subtract for each sweep.
    fs: int or float
        The sampling frequency in kHz. Default is 10 kHz.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    polarity: str
        The expected polarity of the EPSC; negative: '-'; postitive: '+'.
        Default is '-'.
    post_stim: int or float
        Time in ms that marks the end of the sampling window post stimulus.
        Default is 100 ms.
    index: bool
        Determines whether or not to return the peak index in addition to the peak. 


    Returns
    -------
    epsc_peaks: pandas.Series
        The absolute peak of mean baseline subtracted time series data.
    epsc_peak_index: int
        The time at which the peak occurs
    peak_window: pandas.DataFrame
        The window of the time series data where the peak is identified
    '''

    subtracted_data = data - baseline
    start = stim_time * fs
    end = (stim_time + post_stim) * fs
    peak_window = subtracted_data.iloc[start:end]

    if index is True:
        if polarity == '-':
            epsc_peaks = peak_window.min()
            epsc_peaks_index = peak_window.idxmin()
        elif polarity == '+':
            epsc_peaks = peak_window.max()
            epsc_peaks_index = peak_window.idxmax()
        else:
            raise ValueError(
                "polarity must either be + or -"
            )    
        return epsc_peaks, epsc_peaks_index, peak_window

    elif index is False:
        if polarity == '-':
            epsc_peaks = peak_window.min()
        elif polarity == '+':
            epsc_peaks = peak_window.max()
        else:
            raise ValueError(
                "polarity must either be + or -"
            )    
        return epsc_peaks, peak_window


def calculate_latency_jitter(window, stim_time, fs, mean_trace=False):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    window: pandas.Series or pandas.DataFrame
        The time series data window in which the peak is found.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    fs: int or float
        The sampling frequency in kHz.
    mean_trace: bool
        Determines whether or not to return jitter in addition to latency. 

    Returns
    -------
    latency: float or pandas.Series
        The latency to peak from stim onset, in ms
    jitter: float or pandas.Series
        The std of latency to peak (if calculated on multiple traces)
    '''
    peak_time = window.idxmin()
    latency = peak_time - stim_time

    if mean_trace is False:
        jitter = np.std(latency)
        return latency, jitter
    
    elif mean_trace is True:
        return latency



# create dict to hold analysis results for each stim set 
cell_analysis_dict = {}
power_curve_df = pd.DataFrame()

# pick 80%, 2 ms condition to test
# stim_id = 2
# stim_condition = list(sweeps_dict)[stim_id]
# traces = list(sweeps_dict.values())[stim_id]
# mean_trace = traces.mean(axis=1)

for stim_id in range(len(list(sweeps_dict))):
    stim_condition = list(sweeps_dict)[stim_id]
    traces = list(sweeps_dict.values())[stim_id]
    mean_trace = traces.mean(axis=1)
    

    ''' Run the below for each set of stim sweeps '''

    # filter traces
    traces_filtered = filter_traces(traces, fs, lowpass_freq)
    mean_traces_filtered = pd.DataFrame(traces_filtered.mean(axis=1))

    # convert time to ms
    time = np.arange(0, len(traces)/fs, 1/fs)
    mean_traces_filtered.index = time
    traces_filtered.index = time

    # plot mean sweeps to test
    fig = plt.figure()
    plt.plot(mean_traces_filtered)

    # find mean baseline, defined as the last 3s of the sweep
    baseline = calculate_mean_baseline(traces_filtered, fs, baseline_start=3000, baseline_end=6000)
    mean_baseline = calculate_mean_baseline(mean_traces_filtered, fs, baseline_start=3000, baseline_end=6000)

    # find std of baseline
    std_baseline = calculate_std_baseline(traces_filtered, fs, baseline_start=3000, baseline_end=6000)
    mean_std_baseline = calculate_std_baseline(mean_traces_filtered, fs, baseline_start=3000, baseline_end=6000)

    # find current peak
    current_peaks, peak_window = calculate_epsc_peak(traces_filtered, baseline, fs, stim_time, post_stim=100, polarity='-')
    mean_trace_peak, mean_peak_window = calculate_epsc_peak(mean_traces_filtered, mean_baseline, fs, stim_time, post_stim=100, polarity='-')
    current_peaks_mean = current_peaks.mean()

    # find latency to peak, in ms
    latency, jitter = calculate_latency_jitter(peak_window, stim_time, fs, mean_trace=False)
    mean_trace_latency = calculate_latency_jitter(mean_peak_window, stim_time, fs, mean_trace=True)
    latency_mean = latency.mean()

    # collects measurements into cell dict, nested dict for each stim condition
    stim_dict = {}
    stim_dict['Raw Peaks (pA)'] = current_peaks.tolist()
    stim_dict['Mean Raw Peaks (pA)'] = current_peaks_mean
    stim_dict['Mean Trace Peak (pA)'] = mean_trace_peak[0]
    stim_dict['Peak Latencies (ms)'] = latency.tolist()
    stim_dict['Peak Jitter'] = jitter
    stim_dict['Mean Peak Latency (ms)'] = latency_mean
    stim_dict['Mean Trace Latency (ms)'] = mean_trace_latency[0]

    cell_analysis_dict[stim_condition] = stim_dict

    # collects peaks into power_curve_df, column for each stim condition
    stim_peaks = pd.DataFrame(current_peaks)
    stim_peaks.columns = [stim_condition]
    power_curve_df = pd.concat([power_curve_df, stim_peaks], axis=1)


# export analysis values to csv

















