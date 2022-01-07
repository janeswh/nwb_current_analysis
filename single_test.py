import pandas as pd
import matplotlib.pyplot as plt
import pynwb
import os
import numpy as np
from neo.io import IgorIO
from pynwb import NWBHDF5IO
import elephant
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import sem

# gets sweep info for all cells
sweep_info = pd.read_csv('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables/second_dataset_sweep_info.csv', index_col=0)

# reads in NWB file
io = NWBHDF5IO('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data/JH20211103_c3.nwb', 'r', load_namespaces=True)
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
file = 'JH20211103_c3.nwb'
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


def calculate_mean_baseline(data, fs, baseline_start=100, baseline_end=450):
    '''
    Find the mean baseline in a given time series, defined the 100-450 ms window before stimulus onset
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


def calculate_std_baseline(data, fs, baseline_start=100, baseline_end=450):
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


def calculate_latency_jitter(window, epsc_peaks, stim_time, fs, mean_trace=False):
    '''
    Finds the latency of response onset and the jitter of latency.
    Parameters
    ----------
    window: pandas.Series or pandas.DataFrame
        The time series data window in which the peak is found.
    epsc_peaks: float or pandas.Series
        The absolute peak of mean baseline subtracted time series data.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    fs: int or float
        The sampling frequency in kHz.
    mean_trace: bool
        Determines whether or not to return jitter in addition to latency. 

    Returns
    -------
    latency: float or pandas.Series
        The latency of response onset, defined as 5% of epsc_peaks
    jitter: float or pandas.Series
        The std of latency to response onset (if calculated on multiple traces)
    '''
    onset_amp = abs(epsc_peaks * 0.05)
    latency = []
    for sweep in range(len(window.columns)):
        sweep_window = abs(window.iloc[:,sweep])
        onset_idx = np.argmax(sweep_window >= onset_amp[sweep])
        onset_time = sweep_window.index[onset_idx]
        sweep_latency = onset_time - stim_time
        latency.append(sweep_latency)           
    
    if mean_trace is False:
        jitter = np.std(latency)
        return latency, jitter

    elif mean_trace is True:
        return latency
        

def calculate_timetopeak(window, stim_time, fs, mean_trace=False):
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
    time_to_peak = peak_time - stim_time

    return time_to_peak
    

def calculate_responses(baseline_std, peak_mean, timetopeak, threshold=None):
        '''
        Decides on whether there is a response above 2x, 3x above the baseline std,
        or a user-selectable cutoff.
        Parameters
        ----------
        baseline_std: int or float
            The std of the baseline of the mean filtered trace.
        peak_mean: int or float
            The current peak of the mean filtered trace.
        timetopeak: int or float
            The time to current peak. Usually uses the average time to peak, averaged from all sweeps.
        threshold: int, float (optional)
            If supplied, will provide another threshold in addition to the 2x and 3x
            above the baseline std to threshold the response checker.
        Returns
        -------
        esponses: pd.DataFrame(bool)
            A DataFrame with bool for responses above the threshold in the column header.
        '''
        # takes values out of series format to enable boolean comparison
        baseline_std = baseline_std[0]
        peak_mean = peak_mean[0]
        
        response_2x = abs(peak_mean) > baseline_std * 2 and timetopeak < 10
        response_3x = abs(peak_mean) > baseline_std * 3 and timetopeak < 10

        if threshold is None:
            responses = pd.DataFrame({'Response 2x STD': response_2x,
                                           'Response 3x STD': response_3x},
                                           index=range(1))
        else:
            response_threshold = abs(peak_mean) > baseline_std * threshold
            response_string = 'Response {}x STD'.format(threshold)

            responses = pd.DataFrame({'Response 2x STD': response_2x,
                                           'Response 3x STD': response_3x,
                                           response_string: response_threshold},
                                           index=range(1))

        return responses



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
    # mean_trace = traces.mean(axis=1)
    

    ''' Run the below for each set of stim sweeps '''

    # filter traces
    traces_filtered = filter_traces(traces, fs, lowpass_freq)
    mean_trace_filtered = pd.DataFrame(traces_filtered.mean(axis=1))

    # convert time to ms
    time = np.arange(0, len(traces)/fs, 1/fs)
    mean_trace_filtered.index = time
    traces_filtered.index = time

    # plot mean sweeps to test
    fig = plt.figure()
    plt.plot(mean_trace_filtered)

    # find mean baseline, defined as the last 3s of the sweep
    baseline = calculate_mean_baseline(traces_filtered, fs, baseline_start=100, baseline_end=450)
    mean_baseline = calculate_mean_baseline(mean_trace_filtered, fs, baseline_start=100, baseline_end=450)

    # find std of baseline
    std_baseline = calculate_std_baseline(traces_filtered, fs, baseline_start=100, baseline_end=450)
    mean_std_baseline = calculate_std_baseline(mean_trace_filtered, fs, baseline_start=100, baseline_end=450)

    # find current peak
    current_peaks, peak_window = calculate_epsc_peak(traces_filtered, baseline, fs, stim_time, post_stim=100, polarity='-')
    mean_trace_peak, mean_peak_window = calculate_epsc_peak(mean_trace_filtered, mean_baseline, fs, stim_time, post_stim=100, polarity='-')
    current_peaks_mean = current_peaks.mean()


    # find latency to response onset and jitter, in ms
    latency, jitter = calculate_latency_jitter(peak_window, current_peaks, stim_time, fs, mean_trace=False)
    mean_trace_latency = calculate_latency_jitter(mean_peak_window, mean_trace_peak, stim_time, fs, mean_trace=True)
    latency_mean = np.asarray(latency).mean()

    # find time to peak, in ms
    time_to_peak = calculate_timetopeak(peak_window, stim_time, fs, mean_trace=False)
    mean_trace_time_to_peak = calculate_timetopeak(mean_peak_window, stim_time, fs, mean_trace=True)
    time_to_peak_mean = time_to_peak.mean()

    # determines whether the cell is responding, using mean_trace_filtered
    responses = calculate_responses(mean_std_baseline, mean_trace_peak, time_to_peak_mean)

    # collects measurements into cell dict, nested dict for each stim condition
    stim_dict = {}
    stim_dict['Raw Peaks (pA)'] = current_peaks.tolist()
    stim_dict['Mean Raw Peaks (pA)'] = current_peaks_mean
    stim_dict['Mean Trace Peak (pA)'] = mean_trace_peak[0]
    stim_dict['Onset Latencies (ms)'] = latency
    stim_dict['Onset Jitter'] = jitter
    stim_dict['Mean Onset Latency (ms)'] = latency_mean
    stim_dict['Onset SEM'] = sem(latency)
    stim_dict['Mean Trace Onset Latency (ms)'] = mean_trace_latency[0]
    stim_dict['Time to Peaks (ms)'] = time_to_peak.tolist()    
    stim_dict['Mean Time to Peak (ms)'] = time_to_peak_mean
    stim_dict['Time to Peak SEM'] = sem(time_to_peak)
    stim_dict['Mean Trace Time to Peak (ms)'] = mean_trace_time_to_peak[0]


    stim_dict['Response 2x STD'] = responses['Response 2x STD'][0]
    stim_dict['Response 3x STD'] = responses['Response 3x STD'][0]

    cell_analysis_dict[stim_condition] = stim_dict

    # collects peaks into power_curve_df, column for each stim condition
    stim_peaks = pd.DataFrame(current_peaks)
    stim_peaks.columns = [stim_condition]
    power_curve_df = pd.concat([power_curve_df, stim_peaks], axis=1)

# export analysis values to csv



''' The below makes power curve stats and graph '''

# how to convert column names to tuples, then can just pass tuples through to multiindex
tuples = [tuple(condition.split(",")) for condition in power_curve_df.columns]
power_curve_df.columns = pd.MultiIndex.from_tuples(tuples)
power_curve_df.index.names = ['Sweep']
power_curve_df = power_curve_df.T

# get mean response and SEM for plotting
power_curve_stats = pd.DataFrame()
power_curve_stats['Mean Response Amplitude (pA)'] = power_curve_df.mean(axis=1)
power_curve_stats['SEM'] = power_curve_df.sem(axis=1)
power_curve_df = pd.concat([power_curve_df, power_curve_stats], axis=1) # for output to csv

power_curve_stats.reset_index(inplace=True)
power_curve_stats = power_curve_stats.rename(columns={'level_0':'Light Intensity', 'level_1':'Light Duration'})


# plot power curve, with each light pulse duration as a different series
fig = px.line(power_curve_stats, x=power_curve_stats['Light Intensity'], y=power_curve_stats['Mean Response Amplitude (pA)'], 
    color='Light Duration', error_y=power_curve_stats['SEM'], markers=True)
fig['layout']['yaxis']['autorange'] = "reversed"
fig['layout']['xaxis']['autorange'] = "reversed"

fig.show()


''' The below makes response stats plot '''

cell_analysis_df = pd.DataFrame(cell_analysis_dict).T
cell_analysis_df.index = pd.MultiIndex.from_tuples(tuples)



sweep_stats_plt = cell_analysis_df[['Onset Jitter', 'Mean Onset Latency (ms)', 'Mean Trace Onset Latency (ms)',
    'Mean Time to Peak (ms)', 'Mean Trace Time to Peak (ms)']].copy()



# add sem to the various measures


sweep_stats_plt.reset_index(inplace=True)
power_curve_stats = power_curve_stats.rename(columns={'level_0':'Light Intensity', 'level_1':'Light Duration'})


# plot response stats as subplots, with each light pulse duration as a different series
# onset latencies (raw data pts) + SEM
# onset jitter
# mean trace onset latency
# time to peak (raw data pts) + SEM
# mean trace time to peak
# combine with power curve graph?





