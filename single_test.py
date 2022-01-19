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
from plotly.subplots import make_subplots

from scipy.stats import sem

import plotly.io as pio
pio.renderers.default = "browser"

from acq_parameters import *
import pdb



class JaneCell(object):
    def __init__(self, sweep_info, nwbfile):
        self.sweep_info = sweep_info
        self.nwbfile = nwbfile
        self.raw_df = None
        self.sweeps_dict = None

        self.cell_analysis_dict = None
        self.power_curve_df = None
        self.tuples = None
        self.power_curve_stats = None
        self.sweep_analysis_values = None

        self.cell_sweep_info = self.initialize_cell()


    def initialize_cell(self):
        # turn sweep_table into a pandas df
        nwb_df = pd.DataFrame(self.nwbfile.sweep_table[:]) 

        # finds the series indices for VoltageClampSeries data
        b = [series[0] for series in nwb_df['series']]

        vc_indices = []
        for index, item in enumerate(b):
            # print(type(b[index]) == pynwb.icephys.VoltageClampSeries, index)
            if type(b[index]) == pynwb.icephys.VoltageClampSeries:
                vc_indices.append(index)

        # puts the VC series data into another df
        vc_sweeps = nwb_df.iloc[vc_indices]

        # extracts VC series data into df
        raw_series = []

        for i in vc_indices:
            raw_series.append(b[i].data[:])

        # columns are sweeps, rows are time (*fs)
        self.raw_df = pd.DataFrame(raw_series).T
        self.raw_df.columns = vc_sweeps.sweep_number

        # gets sweep info for one cell, drops empty values
        file = 'JH20211103_c3.nwb'
        file_split = file.split('.')
        cell_name = file_split[0]

        cell_sweep_info = pd.DataFrame(self.sweep_info.loc[cell_name]).dropna()

        return cell_sweep_info


    def drop_sweeps(self):
        '''
        If applicable, drops depol sweeps and esc sweeps
        '''

        # pulls the non-escaped sweep numbers out
        if 'Escaped sweeps' in self.cell_sweep_info.index:
            nonesc_range = self.cell_sweep_info.loc['Non-esc sweeps'][0].split(',')
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
        if 'Depol sweeps' in self.cell_sweep_info.index:
            depol_range = self.cell_sweep_info.loc['Depol sweeps'][0].split(',')
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
            self.raw_df = self.raw_df.drop(columns=depol_sweeps, axis=1)

        if 'nonesc_sweeps' in globals():
            self.raw_df = self.raw_df.filter(nonesc_sweeps)

        

    def make_sweeps_dict(self):
        '''
        Create dict with stim name as keys, VC data as values
        '''

        stim_sweep_info = self.cell_sweep_info.filter(like="%", axis=0)
        sweeps_dict = {}

        # define sweep ranges for each stim set present
        for i in range(len(stim_sweep_info.index)):
            stim_name = stim_sweep_info.index[i]
            stim_range = stim_sweep_info.iloc[i].str.split(',')[0]
            stim_sweeps = []

            for j in range(len(stim_range)):
                if '-' in stim_range[j]:
                    r_start = int(stim_range[j].split('-')[0])
                    r_end = int(stim_range[j].split('-')[1])+1
                    all_sweeps = list(range(r_start, r_end))
                    stim_sweeps.extend(all_sweeps)
                else:            
                    stim_sweeps.append(int(stim_range[j][0]))

            stim_sweeps_VC = self.raw_df[self.raw_df.columns.intersection(set(stim_sweeps))]
            sweeps_dict[stim_name] = stim_sweeps_VC

        # drop keys with empty dataframes
        self.sweeps_dict = {k:v for (k,v) in sweeps_dict.items() if not v.empty}

        


    def filter_traces(self, traces):
        '''
        add filtered traces attrbute to data object
        lowpass_freq: frequency in Hz to pass to elephant filter
        '''
        traces_filtered = elephant.signal_processing.butter(traces.T, lowpass_freq=500, fs=FS * 1000)
        traces_filtered = pd.DataFrame(traces_filtered).T

        return traces_filtered


    def calculate_mean_baseline(self, data, baseline_start=100, baseline_end=450):
        '''
        Find the mean baseline in a given time series, defined the 100-450 ms window before stimulus onset
        Parameters
        ----------
        data: pandas.Series or pandas.DataFrame
            The time series data for which you want a baseline.
        baseline_start: int or float
            The time in ms when baseline starts.
        baseline_end: int or float
            The length of the sweep and the time when baseline ends.

        Returns
        -------
        baseline: float or pandas.Series
            The mean baseline over the defined window
        '''
        start = baseline_start * FS
        stop = baseline_end * FS
        window = data.iloc[start:stop]
        baseline = window.mean()

        return baseline


    def calculate_std_baseline(self, data, baseline_start=100, baseline_end=450):
        '''
        Find the mean baseline in a given time series
        Parameters
        ----------
        data: pandas.Series or pandas.DataFrame
            The time series data for which you want a baseline.
        baseline_start: int or float
            The time in ms when baseline starts.
        baseline_end: int or float
            The length of the sweep and the time when baseline ends.

        Returns
        -------
        baseline: float or pandas.Series
            The mean baseline over the defined window
        '''
        start = baseline_start * FS
        stop = baseline_end * FS
        window = data.iloc[start:stop]
        std = window.std()

        return std


    def calculate_epsc_peak(self, data, baseline, post_stim=POST_STIM, polarity='-', index=False):
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
        start = STIM_TIME * FS
        end = (STIM_TIME + post_stim) * FS
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


    def calculate_latency_jitter(self, window, epsc_peaks, mean_trace=False):
        '''
        Finds the latency of response onset and the jitter of latency.
        Parameters
        ----------
        window: pandas.Series or pandas.DataFrame
            The time series data window in which the peak is found.
        epsc_peaks: float or pandas.Series
            The absolute peak of mean baseline subtracted time series data.
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
            sweep_latency = onset_time - STIM_TIME
            latency.append(sweep_latency)           
        
        if mean_trace is False:
            jitter = np.std(latency)
            return latency, jitter

        elif mean_trace is True:
            return latency, None
            

    def calculate_timetopeak(self, window, response_onset, mean_trace=False):
        '''
        Find the mean baseline in a given time series
        Parameters
        ----------
        window: pandas.Series or pandas.DataFrame
            The time series data window in which the peak is found.
        response_onset: int or float
            Time in ms at which response onset occurs (also response latency).
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
        time_to_peak = peak_time - STIM_TIME

        return time_to_peak
        

    def calculate_responses(self, baseline_std, peak_mean, timetopeak, threshold=None):
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


    def calculate_stim_stats(self, stim_id):
        sweeps_dict = self.sweeps_dict

        stim_condition = list(sweeps_dict)[stim_id]
        traces = list(sweeps_dict.values())[stim_id]
        # mean_trace = traces.mean(axis=1)

        ''' Run the below for each set of stim sweeps '''

        # filter traces
        traces_filtered = self.filter_traces(traces)
        mean_trace_filtered = pd.DataFrame(traces_filtered.mean(axis=1))

        # convert time to ms
        time = np.arange(0, len(traces)/FS, 1/FS)
        mean_trace_filtered.index = time
        traces_filtered.index = time

        # plot mean sweeps to test
        fig = plt.figure()
        plt.plot(mean_trace_filtered)

        # find mean baseline, defined as the last 3s of the sweep
        baseline = self.calculate_mean_baseline(traces_filtered, baseline_start=100, baseline_end=450)
        mean_baseline = self.calculate_mean_baseline(mean_trace_filtered, baseline_start=100, baseline_end=450)

        # find std of baseline
        std_baseline = self.calculate_std_baseline(traces_filtered, baseline_start=100, baseline_end=450)
        mean_std_baseline = self.calculate_std_baseline(mean_trace_filtered, baseline_start=100, baseline_end=450)

        # find current peak
        current_peaks, peak_window = self.calculate_epsc_peak(traces_filtered, baseline, post_stim=100, polarity='-')
        mean_trace_peak, mean_peak_window = self.calculate_epsc_peak(mean_trace_filtered, mean_baseline, post_stim=100, polarity='-')
        current_peaks_mean = current_peaks.mean()


        # find latency to response onset and jitter, in ms
        latency, jitter = self.calculate_latency_jitter(peak_window, current_peaks, mean_trace=False)
        mean_trace_latency = self.calculate_latency_jitter(mean_peak_window, mean_trace_peak, mean_trace=True)
        latency_mean = np.asarray(latency).mean()

        # find time to peak, in ms
        time_to_peak = self.calculate_timetopeak(peak_window, latency, mean_trace=False)
        mean_trace_time_to_peak = self.calculate_timetopeak(mean_peak_window, latency, mean_trace=True)
        time_to_peak_mean = time_to_peak.mean()

        # determines whether the cell is responding, using mean_trace_filtered
        responses = self.calculate_responses(mean_std_baseline, mean_trace_peak, time_to_peak_mean)

        # calculates the baseline-subtracted mean trace for plotting purposes
        sub_mean_trace = mean_trace_filtered - mean_baseline

        # collects measurements into cell dict, nested dict for each stim condition
        stim_dict = {
            'Raw Peaks (pA)': current_peaks.tolist(),
            'Mean Raw Peaks (pA)': current_peaks_mean,
            'Mean Trace Peak (pA)': mean_trace_peak[0],
            'Onset Latencies (ms)': latency,
            'Onset Jitter': jitter,
            'Mean Onset Latency (ms)': latency_mean,
            'Onset SEM': sem(latency),
            'Mean Trace Onset Latency (ms)': mean_trace_latency[0][0],
            'Time to Peaks (ms)': time_to_peak.tolist(),
            'Mean Time to Peak (ms)': time_to_peak_mean,
            'Time to Peak SEM': sem(time_to_peak),
            'Mean Trace Time to Peak (ms)': mean_trace_time_to_peak[0],
            'Response 2x STD': responses['Response 2x STD'][0],
            'Response 3x STD': responses['Response 3x STD'][0]
        }

        return stim_condition, sub_mean_trace, current_peaks, stim_dict


    def make_cell_analysis_dict(self):

        # create dict to hold analysis results for each stim set 
        cell_analysis_dict = {}
        power_curve_df = pd.DataFrame()
        all_mean_traces = pd.DataFrame()

        for stim_id in range(len(list(self.sweeps_dict))):
            stim_condition, sub_mean_trace, current_peaks, stim_dict = self.calculate_stim_stats(stim_id)
            cell_analysis_dict[stim_condition] = stim_dict

            # collects peaks into power_curve_df, column for each stim condition
            stim_peaks = pd.DataFrame(current_peaks)
            stim_peaks.columns = [stim_condition]
            power_curve_df = pd.concat([power_curve_df, stim_peaks], axis=1)

            # collects mean traces into all_mean_traces, column for each stim condition
            mean_trace = pd.DataFrame(sub_mean_trace)
            mean_trace.columns = [stim_condition]
            all_mean_traces = pd.concat([all_mean_traces, mean_trace], axis=1)

        self.cell_analysis_dict = cell_analysis_dict
        self.power_curve_df = power_curve_df
        self.all_mean_traces = all_mean_traces

    
    # export analysis values to csv


    def make_power_curve_stats_df(self):

        ''' 
        The below makes power curve stats table used to plot power curve
        '''
        power_curve_df = self.power_curve_df

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

        self.tuples = tuples
        self.power_curve_stats = power_curve_stats


    def make_stats_df(self):

        ''' 
        The below makes response stats tables used to plot graphs
        '''
        cell_analysis_df = pd.DataFrame(self.cell_analysis_dict).T
        cell_analysis_df.index = pd.MultiIndex.from_tuples(self.tuples)


        sweep_analysis_values =  cell_analysis_df[['Onset Latencies (ms)', 'Time to Peaks (ms)']].copy()
        sweep_analysis_values = sweep_analysis_values.explode(['Onset Latencies (ms)', 'Time to Peaks (ms)'])

        sweep_analysis_values.reset_index(inplace=True)
        sweep_analysis_values = sweep_analysis_values.rename(columns={'level_0':'Light Intensity', 'level_1':'Light Duration'})

        cell_analysis_df.reset_index(inplace=True)
        cell_analysis_df = cell_analysis_df.rename(columns={'level_0':'Light Intensity', 'level_1':'Light Duration'})
        
        self.sweep_analysis_values = sweep_analysis_values
        self.cell_analysis_df = cell_analysis_df
        

    def graph_curve_stats(self):
        '''
        do a loop through available durations and intensities instead of hard
        coding. maybe need MultiIndex after all?? Put power curve + all stats 
        measurements in subplots
        '''

        power_curve_stats = self.power_curve_stats
        sweep_analysis_values = self.sweep_analysis_values
        cell_analysis_df = self.cell_analysis_df

        intensities = sweep_analysis_values['Light Intensity'].unique()
        durations = sweep_analysis_values['Light Duration'].unique()


        color = px.colors.qualitative.Plotly
        fig = make_subplots(
            rows=3, cols=2,
            x_title = "Light Intensity (%)"
            )



        for count, duration in enumerate(durations):

            error = power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['SEM']].squeeze()
            
            # power curve
            fig.add_trace(go.Scatter(
                x=power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['Mean Response Amplitude (pA)']].squeeze(), 
                name=duration,
                error_y=dict(
                    type='data',
                    array=error.values,
                    visible=True),     
                mode='lines+markers',
                line=dict(color=color[count])
                ),
                row=1, col=1
            )
            
            # onset latency
            fig.add_trace(go.Box(
                x=sweep_analysis_values.loc[sweep_analysis_values['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=sweep_analysis_values.loc[sweep_analysis_values['Light Duration']==duration,['Onset Latencies (ms)']].squeeze(),
                name=duration,
                line=dict(color=color[count]),
                showlegend=False,
                ),
                row=1, col=2
            )
            
            # onset jitter
            fig.add_trace(go.Bar(
                x=power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=cell_analysis_df.loc[cell_analysis_df['Light Duration']==duration, ['Onset Jitter']].squeeze(),
                name=duration,
                marker_color=color[count],
                showlegend=False,
                ),
                row=2, col=1
            )
            
            
            # mean trace onset latency
            fig.add_trace(go.Bar(
                x=power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=cell_analysis_df.loc[cell_analysis_df['Light Duration']==duration, ['Mean Trace Onset Latency (ms)']].squeeze(),
                name=duration,
                marker_color=color[count],
                showlegend=False,
                ),
                row=2, col=2
            )


            # time to peak
            fig.add_trace(go.Box(
                x=sweep_analysis_values.loc[sweep_analysis_values['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=sweep_analysis_values.loc[sweep_analysis_values['Light Duration']==duration,['Time to Peaks (ms)']].squeeze(),
                name=duration,
                line=dict(color=color[count]),
                showlegend=False,
                ),
                row=3, col=1
            )

            # mean trace time to peak
            fig.add_trace(go.Bar(
                x=power_curve_stats.loc[power_curve_stats['Light Duration']==duration,['Light Intensity']].squeeze(),
                y=cell_analysis_df.loc[cell_analysis_df['Light Duration']==duration, ['Mean Trace Time to Peak (ms)']].squeeze(),
                name=duration,
                marker_color=color[count],
                showlegend=False,
                ),
                row=3, col=2
            )

            # Update xaxis properties
            fig.update_xaxes(autorange='reversed')


            # Update yaxis properties
            fig.update_yaxes(title_text="Response Amplitude (pA)", row=1, col=1, autorange='reversed')
            fig.update_yaxes(title_text="Onset Latency (ms)", row=1, col=2)
            fig.update_yaxes(title_text="Onset Jitter", row=2, col=1)
            fig.update_yaxes(title_text="Mean Trace Onset Latency (ms)", row=2, col=2)
            fig.update_yaxes(title_text="Time to Peak (ms)", row=3, col=1)
            fig.update_yaxes(title_text="Mean Trace Time to Peak (ms)", row=3, col=2)

        fig.update_layout(
            #yaxis_title='Onset Latency (ms)',
            boxmode='group' # group together boxes of the different traces for each value of x
        )

        fig.show()


    def make_mean_traces_df(self):

        ''' 
        The below makes mean traces df used to plot graphs
        '''
        
        mean_trace_df = self.all_mean_traces.T
        mean_trace_df.index = pd.MultiIndex.from_tuples(self.tuples)

        mean_trace_df.reset_index(inplace=True)
        mean_trace_df = mean_trace_df.rename(columns={'level_0':'Light Intensity', 'level_1':'Light Duration'})

        self.mean_trace_df = mean_trace_df


    def graph_response_trace(self):
        '''
        Plots the baseline-subtracted mean trace for each stimulus condition around the response time, 
        one subplot for each duration, if applicable
        '''

        # intensities and durations, and color might need to become self variables

        sweep_analysis_values = self.sweep_analysis_values
        intensities = sweep_analysis_values['Light Intensity'].unique()
        durations = sweep_analysis_values['Light Duration'].unique()

        color = px.colors.qualitative.Plotly

        

        stim_columns = self.mean_trace_df.loc[:,['Light Intensity', 'Light Duration']]
        traces_to_plot = self.mean_trace_df.loc[:,500.00:700.00] # only plots first 400-1000 ms
        traces_to_plot_combined = pd.concat([stim_columns, traces_to_plot], axis=1)


        fig = make_subplots(
            #rows=len(intensities), cols=1,
            rows=1, cols=len(intensities),
            subplot_titles=(intensities + " Light Intensity"),
            shared_yaxes=True,
            x_title = "Time (ms)",
            y_title = "Amplitude (pA)"
            )

        

        for intensity_count, intensity in enumerate(intensities):
            for duration_count, duration in enumerate(durations):

                # plot sweeps from all intensities of one duration
                fig.add_trace(go.Scatter(
                    x=traces_to_plot.columns,
                    y=traces_to_plot_combined.loc[(traces_to_plot_combined['Light Intensity']==intensity) 
                        & (traces_to_plot_combined['Light Duration']==duration), 500.00::].squeeze(), 
                    name=duration,
                    mode='lines',
                    line=dict(color=color[duration_count]),
                    legendgroup=duration,
                    ),
                    #row=intensity_count+1, col=1
                    row=1, col=intensity_count+1
                )

        # below is code from stack overflow to hide duplicate legends
        names = set()
        fig.for_each_trace(
            lambda trace:
                trace.update(showlegend=False)
                if (trace.name in names) else names.add(trace.name))
        
        fig.update_layout(legend_title_text='Light Duration')

        pdb.set_trace()

        

        fig.show()



        
        # fig2 = px.line(trace_to_plot, x=trace_to_plot.index, y=trace_to_plot.columns)
        # fig2.show()


    

