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


def run_single(dataset, csvfile, nwbfile_name):
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
    cell.make_drug_sweeps_dict()
    cell.make_sweeps_dict()
    cell.extract_drug_sweeps()

    # 4 runs stats on sweeps and creates a dict for each stim condition
    cell.make_cell_analysis_dict()

    # 5 calculates power curve for plotting
    cell.make_power_curve_stats_df()

    # 6 calculates response stats for plotting
    cell.make_stats_df()

    # 7 plots mean traces for ctrl vs. NBQX wash-in (from the same cell)
    type = "drug"
    traces = cell.make_mean_traces_df()
    drug_trace = cell.extract_drug_sweeps()
    axes, noaxes = plot_example_traces(traces, drug_trace, type)
    save_example_traces_figs(axes, noaxes, type)

    # cell.graph_response_trace()

    # # 8 makes plots for power curve and response stats if cell responds
    # if response == True:

    #     summary_plots = cell.graph_curve_stats()
    #     cell.export_stats_csv()
    # else:
    #     print("Cell doesn't have response, no response stats plotted")

    # # 9 saves combined plots as html file, exports stats as csv
    # cell.output_html_plots()


if __name__ == "__main__":
    dataset = "non-injected"
    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
        dataset,
        csvfile_name,
    )
    nwbfile_name = "JH20211130_c1.nwb"

    run_single(dataset, csvfile, nwbfile_name)

    print("Analysis for {} done".format(nwbfile_name))
