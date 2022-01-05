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

baseline_start = 3000
baseline_end = 6000

stim_conditions = list(sweeps_dict)
genotype = cell_sweep_info.loc['Genotype']

''' Functions '''

def filter_traces(traces, fs, lowpass_freq):
    '''
    add filtered traces attrbute to data object
    lowpass_freq: frequency in Hz to pass to elephant filter
    '''
    traces_filtered = elephant.signal_processing.butter(traces.T, lowpass_freq=lowpass_freq, fs=self.fs * 1000)
    traces_filtered = pd.DataFrame(traces_filtered).T

    return traces_filtered


''' Run the below for each set of stim sweeps, using mean filtered traces below '''

# pick 80%, 2 ms condition to test
traces = list(sweeps_dict.values())[2]
mean_traces = traces.mean(axis=1)

time = np.arange(0, len(traces)/fs, 1/fs)

# filter traces

mean_traces_filtered = filter_traces(traces, fs, lowpass_freq)
traces_filtered = elephant.signal_processing.butter(traces.T, lowpass_freq=lowpass_freq, fs=fs * 1000)
traces_filtered = pd.DataFrame(traces_filtered).T
mean_traces_filtered = pd.DataFrame(traces_filtered.mean(axis=1))
mean_traces_filtered.index = time

# plot sweeps to test
fig = plt.figure()
plt.plot(mean_traces_filtered)

# find mean baseline, defined as the last 3s of the sweep
start = baseline_start * fs
stop = baseline_end * fs
window = mean_traces_filtered.iloc[start:stop]
baseline = window.mean()

# find std of baseline
start = baseline_start * fs
stop = baseline_end * fs
window = mean_traces_filtered.iloc[start:stop]
std = window.std()

# find current peak
subtracted_data = mean_traces_filtered - baseline
post_stim = 100 # 100 ms window after stimulus for finding peak
start = stim_time * fs
end = (stim_time + post_stim) * fs
peak_window = subtracted_data.iloc[start:end]
polarity = "-"

if polarity == '-':
    epsc_peaks = peak_window.min()
elif polarity == '+':
    epsc_peaks = peak_window.max()
else:
    raise ValueError(
        "polarity must either be + or -"
    )    

# find latency to peak (one latency from mean trace), in ms
peak_time = peak_window.idxmin()
latency = peak_time - stim_time

# find jitter by calculating all the latencies in the set of stim sweeps


















# import clamp_ephys
# import pandas as pd
# import os
# import scipy
# import matplotlib
# import matplotlib.pyplot as plt

# '''####################### SET THE PROPER PATH YOU WANT ########################### '''
# paths = clamp_ephys.workflows.file_structure('local', 'MMZ_STC_dataset')
# tables = paths.tables
# figures = paths.figures











# data_raw = IgorIO(filename='/home/jhuang/Documents/phd_projects/Injected_GC_data_060820/VC_pairs/data/p2/JH200303_c1_light100.ibw')
# data_neo = data_raw.read_block()
# data_neo_array = data_neo.segments[0].analogsignals[0]
# data_df = pd.DataFrame(data_neo_array.as_array().squeeze())

















