import pandas as pd
import pynwb
import os
import numpy as np
from neo.io import IgorIO
from pynwb import NWBHDF5IO

# reads in NWB file
io = NWBHDF5IO('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data/JH20211005_c3.nwb', 'r', load_namespaces=True)
nwbfile = io.read()

# gets the IC timeseries data from the first sweep
nwbfile.sweep_table.series[0].data[:] 

# turn sweep_table into a pandas df
nwb_df = pd.DataFrame(nwbfile.sweep_table[:]) 

# defines list of sweeps to pull from sweep table, then extracts into subset
select_sweeps = list(range(11, 21)) 
subset_nwb_df = nwb_df[nwb_df['sweep_number'].isin(select_sweeps)] 

# finds the series indices for VoltageClampSeries data
b = [series[0] for series in subset_nwb_df['series']]

vc_indices = []
for index, item in enumerate(b):
    print(type(b[index]) == pynwb.icephys.VoltageClampSeries, index)
    if type(b[index]) == pynwb.icephys.VoltageClampSeries:
        vc_indices.append(index)

# puts the VC series data into another df
vc_sweeps = subset_nwb_df.iloc[vc_indices]

# extracts VC series data into df
raw_series = []

for i in vc_indices:
    raw_series.append(b[i].data[:])

# columns are sweeps, rows are time (*fs)
raw_df = pd.DataFrame(raw_series).T
raw_df.columns = vc_sweeps.sweep_number





data_raw = IgorIO(filename='/home/jhuang/Documents/phd_projects/Injected_GC_data_060820/VC_pairs/data/p2/JH200303_c1_light100.ibw')
data_neo = data_raw.read_block()
data_neo_array = data_neo.segments[0].analogsignals[0]
data_df = pd.DataFrame(data_neo_array.as_array().squeeze())


















subset_nwb_df.iloc[0][0][0].data[:] # the VoltageClampSeries from one sweep


subset_nwb_df['series'].values[:][0][0].data[:]


type(subset_nwb_df['series'][50][0])








sweep159 = nwb_df.loc[nwb_df['sweep_number']==159] # select all the series from sweep 159
sweep159.loc[466][0]    # trying to isolate timeseries data



nwb_df.loc[nwb_df['series'].str.contains('VoltageClampSeries')]