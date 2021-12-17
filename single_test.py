import pandas as pd
import pynwb
import os
import numpy as np
from neo.io import IgorIO
from pynwb import NWBHDF5IO

# gets sweep info for all cells
sweep_info = pd.read_csv('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables/second_dataset_sweep_info.csv', index_col=0)


# gets sweep info for one cell, drops empty values
file = 'JH20211006_c1'
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

stim_sweep_info = cell_sweep_info.filter(like="%", axis=0)

# create dict with stim name as keys, VC data as values
sweeps_dict = {}

# define sweep ranges for each stim set present
for i in range(len(stim_sweep_info.index)):
    stim_name = stim_sweep_info.index[i]
    stim_range = stim_sweep_info.iloc[i].str.split(',')[0]
    stim_sweeps = []
    if len(stim_range) == 1:
        if '-' not in stim_range[0][0]:
            stim_sweeps.append(int(stim_range[0][0]))
        else:
            r_start = int(stim_range[0][0].split('-')[0])
            r_end = int(stim_range[0][0].split('-')[1])+1
            all_sweeps = list(range(r_start, r_end))
            stim_sweeps.extend(all_sweeps)
    else:
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



# reads in NWB file
io = NWBHDF5IO('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data/JH20211006_c1.nwb', 'r', load_namespaces=True)
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


# if applicable, drops depol sweeps and  esc sweeps
if 'depol_sweeps' in globals():
    raw_df = raw_df.drop(columns=depol_sweeps, axis=1)

if 'nonesc_sweeps' in globals():
    raw_df = raw_df.filter(nonesc_sweeps)





























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

















