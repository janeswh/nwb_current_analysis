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
from file_settings import FileSettings

from scipy.stats import sem

import plotly.io as pio

pio.renderers.default = "browser"

from acq_parameters import *
import pdb
from single_test import JaneCell
from plotting import *


def get_single_cell(dataset, csvfile, nwbfile_name):
    """
    Initiates the cell object for extracting traces later
    """
    nwbfile = os.path.join(
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data",
        dataset,
        nwbfile_name,
    )

    # gets sweep info for all cells
    sweep_info = pd.read_csv(csvfile, index_col=0)

    # reads in NWB file
    io = NWBHDF5IO(nwbfile, "r", load_namespaces=True)
    nwbfile = io.read()

    # 0 initializes JaneCell class
    cell = JaneCell(dataset, sweep_info, nwbfile, nwbfile_name)

    return cell


# def get_cell_sweeps_dict(cell, spikes=False):
#     """
#     Gets the dict of all sweeps for a cell object
#     """
#     # 2 drops depolarized and esc AP sweeps from VC data if applicable
#     cell.drop_sweeps()

#     # 3 makes a dict for each cell, with stim condition as keys and all sweeps per stimulus as values
#     if spikes is False:
#         cell.make_sweeps_dict()
#         return cell.sweeps_dict

#     else:
#         cell.make_spikes_dict()
#         return cell.ic_sweeps_dict


def get_annotation_values(cell, selected_condition, sweep_number):
    """
    Gets the response onset, amplitude peak, and time to peak values for
    annotating example trace.
    """

    onset_time = cell.cell_analysis_dict[selected_condition][
        "Onset Times (ms)"
    ][sweep_number]

    onset_latency = cell.cell_analysis_dict[selected_condition][
        "Onset Latencies (ms)"
    ][sweep_number]

    amplitude_peak = cell.cell_analysis_dict[selected_condition][
        "Raw Peaks (pA)"
    ][sweep_number]

    time_topeak = cell.cell_analysis_dict[selected_condition][
        "Time to Peaks (ms)"
    ][sweep_number]

    return [onset_time, onset_latency, amplitude_peak, time_topeak]


def get_single_cell_traces(
    cell, traces_type="mean", sweep_number=None, annotate=False
):
    """
    Gets the normal light stim mean traces for a cell object
    """

    # 2 drops depolarized and esc AP sweeps from VC data if applicable
    cell.drop_sweeps()

    # 3 makes a dict for each cell, with stim condition as keys and all sweeps per stimulus as values
    if traces_type == "spikes":
        cell.make_spikes_dict()
    elif traces_type == "oscillations":
        cell.make_resting_dict()
    else:
        cell.make_sweeps_dict()

        # 4 runs stats on sweeps and creates a dict for each stim condition
        cell.make_cell_analysis_dict()

        # 5 calculates power curve for plotting
        cell.make_power_curve_stats_df()

        # 6 calculates response stats for plotting
        cell.make_stats_df()

    if cell.cell_name == "JH20210923_c2":
        selected_condition = "50%, 1 ms"
    else:
        selected_condition = ",".join(FileSettings.SELECTED_CONDITION)

    if annotate is True:
        annotation_values = get_annotation_values(
            cell, selected_condition, sweep_number
        )
    else:
        annotation_values = None

    if traces_type == "mean":
        cell.make_mean_traces_df()
        traces = cell.mean_trace_df

    elif traces_type == "spikes":
        # # pulls out spikes sweep
        spike_sweep = cell.extract_FI_sweep(sweep_number)
        traces = spike_sweep

    elif traces_type == "single":
        vc_sweep = cell.filtered_traces_dict[selected_condition][sweep_number]
        # vc_sweep = cell.extract_VC_sweep(selected_condition, sweep_number)
        traces = vc_sweep

    return traces, annotation_values


def get_single_drug_traces(cell):
    """
    Gets the drug traces for a cell object
    """

    cell.make_drug_sweeps_dict()
    drug_trace = cell.extract_drug_sweeps()

    return drug_trace


def make_example_traces(dataset, csvfile, genotype, main_plot_files):
    """
    Makes plotting traces without inset plots.
    """
    if dataset == "dox_5dpi":
        small_scalebar = True
    else:
        small_scalebar = False

    main_type_names = [
        "{} cell 1".format(genotype),
        "{} cell 2".format(genotype),
    ]
    ephys_traces_plotted = pd.DataFrame()
    # gets the ephys traces for main plot cells
    main_ephys_traces = []
    for file in main_plot_files:
        cell = get_single_cell(dataset, csvfile, file)
        traces, annotation_values = get_single_cell_traces(cell)
        main_ephys_traces.append(traces)

    # makes the plotting traces for main plot
    main_plot_traces = []
    for count, trace in enumerate(main_ephys_traces):
        plot_trace = make_one_plot_trace(
            main_plot_files[count], trace, main_type_names[count]
        )
        main_plot_traces.append(plot_trace)
        plotted_trace = pd.DataFrame(
            {main_type_names[count]: plot_trace["y"]}, index=plot_trace["x"]
        )
        ephys_traces_plotted = pd.concat(
            [ephys_traces_plotted, plotted_trace], axis=1
        )
    fig = plot_example_traces(
        genotype, main_plot_traces, small_scalebar=small_scalebar
    )
    new_save_example_traces_figs(fig, ephys_traces_plotted, dataset)
    # save_fig_to_png(
    #     fig,
    #     legend=True,
    #     rows=1,
    #     cols=2,
    #     png_filename="3dpi_MMZ_example_traces.png",
    # )


def make_inset_plot(
    dataset, csvfile, genotype, main_plot_files, inset_plot_file
):
    """
    Gets the ephys traces and passes them to make plotting traces, then make
    inset plot.
    """

    main_type_names = [
        "{} cell 1".format(genotype),
        "{} cell 2".format(genotype),
    ]

    # gets the ephys traces for main plot cells
    main_ephys_traces = []
    for file in main_plot_files:
        cell = get_single_cell(dataset, csvfile, file)
        traces, annotation_values = get_single_cell_traces(cell)
        main_ephys_traces.append(traces)

    # gets ephys traces for inest plot cell
    inset_cell = get_single_cell(dataset, csvfile, inset_plot_file)
    inset_traces, annotation_values = get_single_cell_traces(inset_cell)
    drug_trace = get_single_drug_traces(inset_cell)

    # makes the plotting traces for main plot
    main_plot_traces = []
    for count, trace in enumerate(main_ephys_traces):
        plot_trace = make_one_plot_trace(
            main_plot_files[count], trace, main_type_names[count]
        )

        main_plot_traces.append(plot_trace)

    # makes the plotting traces for inset plot
    inset_ctrl_trace = make_one_plot_trace(
        inset_plot_file, inset_traces, "Control", inset=True
    )
    inset_drug_trace = make_one_plot_trace(
        inset_plot_file, drug_trace, "NBQX", inset=True
    )

    # puts everything in main plot + inset
    axes, noaxes = make_inset_plot_fig(
        genotype,
        main_plot_traces[0],
        main_plot_traces[1],
        inset_ctrl_trace,
        inset_drug_trace,
    )

    # saves figs
    save_example_traces_figs(axes, noaxes, genotype)

    print("Finished saving inset plots")


def make_annotated_trace(dataset, csvfile, genotype, file_name, sweep_number):
    """
    Plots a single VC trace to demonstrate onset latency, peak amplitude, and
    time to peak.
    """

    cell = get_single_cell(dataset, csvfile, file_name)
    traces, annotation_values = get_single_cell_traces(
        cell, "single", sweep_number, annotate=True
    )
    axes, noaxes = plot_annotated_trace(traces, annotation_values, genotype)
    save_annotated_figs(axes, noaxes, cell, genotype)

    print("Finished saving annotated trace plots")


def make_spike_traces(dataset, csvfile, genotype, file_name, sweep_number):
    """
    Plots a single IC trace to demonstrate STC spike shapes, then also plots
    zoomed in version of the first few spikes.
    """
    cell = get_single_cell(dataset, csvfile, file_name)
    trace, annotation_values = get_single_cell_traces(
        cell, traces_type="spikes", sweep_number=sweep_number
    )
    axes, noaxes = plot_spike_sweeps(genotype, trace)
    pdb.set_trace()
    save_spike_figs(axes, noaxes, cell, genotype)

    print("Finished saving spike plots")


def make_oscillation_trace(
    dataset, csvfile, genotype, file_name, sweep_number
):
    """
    Plots a single IC trace to demonstrate STC oscillation, uses ibw
    """
    path_to_file = os.path.join(FileSettings.DATA_FOLDER, "extra", file_name)
    data_raw = IgorIO(filename=path_to_file)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array().squeeze())

    # filter traces
    traces_filtered = elephant.signal_processing.butter(
        data_df.T, lowpass_freq=500, fs=FS * 1000
    )

    traces_filtered = pd.DataFrame(traces_filtered).T

    time = np.arange(0, len(data_df) / FS, 1 / FS)
    traces_filtered.index = time
    sweep_to_plot = traces_filtered[sweep_number]
    fig = plot_oscillation_sweep(sweep_to_plot)
    new_save_example_traces_figs(fig, sweep_to_plot, "oscillation")

    pdb.set_trace()

    print("Finished saving oscillation plot")


def make_power_curves(dataset, csvfile, genotype, file_name):
    """
    Plots example traces and power curve amplitudes for one cell
    """
    cell = get_single_cell(dataset, csvfile, file_name)
    # 2 drops depolarized and esc AP sweeps from VC data if applicable
    cell.drop_sweeps()

    # 3 makes a dict for each cell, with stim condition as keys and all sweeps per stimulus as values
    cell.make_sweeps_dict()

    # 4 runs stats on sweeps and creates a dict for each stim condition
    cell.make_cell_analysis_dict()

    # 5 calculates power curve for plotting
    cell.make_power_curve_stats_df()

    # 6 calculates response stats for plotting
    cell.make_stats_df()

    # 7 plots mean traces
    cell.make_mean_traces_df()
    power_curve_traces = plot_power_curve_traces(
        cell.mean_trace_df, cell.sweep_analysis_values
    )

    # change only the Time (ms) annotation to have more room
    power_curve_traces.layout.annotations[5].yshift = -50
    save_power_curve_traces(genotype, cell.cell_name, power_curve_traces)

    power_curve_fig = graph_power_curve(
        cell.power_curve_stats, cell.sweep_analysis_values
    )

    save_power_curve(genotype, cell.cell_name, power_curve_fig)


if __name__ == "__main__":
    dataset = "non-injected"
    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
        dataset,
        csvfile_name,
    )

    # # inset plot for Gg8, list big response cell first
    # main_plot_files = ["JH20210923_c2.nwb", "JH20210922_c1.nwb"]
    # inset_plot_file = "JH20211130_c1.nwb"
    # make_inset_plot(dataset, csvfile, "Gg8", main_plot_files, inset_plot_file)

    # # inset plot for OMP
    # main_plot_files = ["JH20211005_c3.nwb", "JH20211029_c1.nwb"]
    # inset_plot_file = "JH20211103_c3.nwb"
    # make_inset_plot(dataset, csvfile, "OMP", main_plot_files, inset_plot_file)

    # # plot single VC trace to show onset latency, pick sweep 131
    # make_annotated_trace(dataset, csvfile, "Gg8", "JH20210923_c2.nwb", 0)

    # # plot one IC trace to show STC spikes, JH20211130_c1 sweep 4 (Gg8)
    # make_spike_traces(dataset, csvfile, "Gg8", "JH20211130_c1.nwb", 4)

    # # plot one IC trace to show STC oscillations, JH20211007_c1 sweep 1 (OMP)
    # make_oscillation_trace(
    #     dataset, csvfile, "Gg8", "JH202111007_c1_oscillation.ibw", 1
    # )

    # # plot example response traces and power curve amplitudes for one OMP cell
    # # JH20211103_c3
    # make_power_curves(dataset, csvfile, "OMP", "JH20211103_c3.nwb")

    # # example traces for 3dpi MMZ

    # dataset = "3dpi"
    # csvfile_name = "{}_sweep_info.csv".format(dataset)
    # csvfile = os.path.join(
    #     "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
    #     dataset,
    #     csvfile_name,
    # )
    # main_plot_files = ["JH20211202_c1.nwb", "JH20211202_c2.nwb"]
    # make_example_traces(dataset, csvfile, "3 dpi MMZ", main_plot_files)

    # example traces for dox 5dpi MMZ
    dataset = "dox_5dpi"
    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
        dataset,
        csvfile_name,
    )
    main_plot_files = ["JH20220111_c4.nwb"]
    make_example_traces(dataset, csvfile, "Dox 5 dpi", main_plot_files)
