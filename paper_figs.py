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


def plot_single_cell(dataset, csvfile, nwbfile_name, type):
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

    # 1 checks whether cell has a response before proceeding
    response = cell.check_response()

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
    traces = cell.make_mean_traces_df()

    # 7 plots mean traces for ctrl vs. NBQX wash-in (from the same cell)
    if type == "drug":
        cell.make_drug_sweeps_dict()
        drug_trace = cell.extract_drug_sweeps()
        axes, noaxes = plot_example_traces(traces, drug_trace, type)
        save_example_traces_figs(axes, noaxes, type, cell.cell_name)

    if type == "genotype":
        return traces


def plot_two_cells(dataset, csvfile, nwbfile_names, type):
    traces = []
    for file in nwbfile_names:
        trace = plot_single_cell(dataset, csvfile, file, type)
        traces.append(trace)

    if nwbfile_names[1] == "JH20210923_c2.nwb":
        exception = True
    else:
        exception = False

    axes, noaxes = plot_example_traces(traces[0], traces[1], type, exception)
    save_example_traces_figs(axes, noaxes, type)


if __name__ == "__main__":
    dataset = "non-injected"
    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
        dataset,
        csvfile_name,
    )

    # plot ctrl vs NBQX traces for one OMP cell
    nwbfile_name = "JH20211103_c3.nwb"
    plot_single_cell(dataset, csvfile, nwbfile_name, "drug")
    print("Analysis for {} done".format(nwbfile_name))

    # plot ctrl vs NBQX traces for one Gg8 cell
    nwbfile_name = "JH20211130_c1.nwb"
    plot_single_cell(dataset, csvfile, nwbfile_name, "drug")
    print("Analysis for {} done".format(nwbfile_name))

    # # plot OMP vs Gg8 traces for two cells - close to median values
    # nwbfile_names = ["JH20211029_c1.nwb", "JH20210922_c1.nwb"]
    # plot_two_cells(dataset, csvfile, nwbfile_names, "genotype")

    # plot OMP vs Gg8 traces for two cells - ideal cells (big responses)
    # actually this is complicated because don't have 100%, 1 ms, only 50% at
    # 1 ms so would need to index differently/make exception
    nwbfile_names = ["JH20211103_c3.nwb", "JH20210923_c2.nwb"]
    plot_two_cells(dataset, csvfile, nwbfile_names, "genotype")

    print("plotting done")

