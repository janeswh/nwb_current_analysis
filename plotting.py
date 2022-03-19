import pandas as pd
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import sem
from collections import defaultdict
import plotly.io as pio
from file_settings import FileSettings

pio.renderers.default = "browser"
import pdb


def plot_averages(dataset, genotype, threshold, averages_df):

    # plotting
    durations = averages_df["Light Duration"].unique()
    intensities = averages_df["Light Intensity"].unique()
    color_dict = {
        " 2 ms": "#7D1935",
        " 1 ms": "#B42B51",
        " 0.25 ms": "#E63E6D",
        " 0.01 ms": "#F892B9",
    }
    summary_stats_fig = make_subplots(
        rows=3, cols=2, x_title="Light Intensity (%)"
    )
    # pdb.set_trace()
    x_intensity = averages_df["Light Intensity"].unique()

    for count, duration in enumerate(durations):
        # mean trace peak amplitude
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Mean Trace Peak (pA)"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=1,
            col=1,
        )

        # Mean Onset Latency (ms)
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Mean Onset Latency (ms)"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=1,
            col=2,
        )

        # onset jitter
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Onset Jitter"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=2,
            col=1,
        )

        # mean trace onset latency
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Mean Trace Onset Latency (ms)"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=2,
            col=2,
        )

        # mean time to peak
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Mean Time to Peak (ms)"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=3,
            col=1,
        )

        # mean trace time to peak
        summary_stats_fig.add_trace(
            go.Box(
                x=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Light Intensity"],
                ].squeeze(),
                y=averages_df.loc[
                    averages_df["Light Duration"] == duration,
                    ["Mean Trace Time to Peak (ms)"],
                ].squeeze(),
                name=duration,
                line=dict(color=color_dict[duration]),
                legendgroup=duration,
            ),
            row=3,
            col=2,
        )

        # Update xaxis properties
        # summary_stats_fig.update_xaxes(autorange="reversed")
        # this defines the intensities order for x-axes
        summary_stats_fig.update_xaxes(
            categoryorder="array", categoryarray=np.flip(intensities)
        )

        # Update yaxis properties
        summary_stats_fig.update_yaxes(
            title_text="Mean Response Amplitude (pA)",
            row=1,
            col=1,
            autorange="reversed",
        )
        summary_stats_fig.update_yaxes(
            title_text="Mean Onset Latency (ms)", row=1, col=2
        )
        summary_stats_fig.update_yaxes(
            title_text="Mean Onset Jitter", row=2, col=1
        )
        summary_stats_fig.update_yaxes(
            title_text="Mean Trace Onset Latency (ms)", row=2, col=2
        )
        summary_stats_fig.update_yaxes(
            title_text="Mean Time to Peak (ms)", row=3, col=1
        )
        summary_stats_fig.update_yaxes(
            title_text="Mean Trace Time to Peak (ms)", row=3, col=2
        )

    summary_stats_fig.update_layout(
        # yaxis_title='Onset Latency (ms)',
        boxmode="group"  # group together boxes of the different traces for each value of x
    )

    # below is code from stack overflow to hide duplicate legends
    names = set()
    summary_stats_fig.for_each_trace(
        lambda trace: trace.update(showlegend=False)
        if (trace.name in names)
        else names.add(trace.name)
    )

    summary_stats_fig.update_layout(
        legend_title_text="Light Duration",
        title_text=(
            dataset
            + " "
            + genotype
            + " summary values, "
            + str(threshold)
            + " mean onset latency threshold"
        ),
        title_x=0.5,
    )

    # summary_stats_fig.show()

    return summary_stats_fig


def save_summary_stats_fig(genotype, threshold, fig_folder, fig):

    html_filename = "{}_{}_threshold_summary_avgs.html".format(
        genotype, threshold
    )
    path = os.path.join(fig_folder, html_filename)

    fig.write_html(path, full_html=False, include_plotlyjs="cdn")


def plot_selected_averages(threshold, selected_avgs_dict):

    # do this for each threshold
    selected_avgs = selected_avgs_dict[threshold]

    genotype_color = {"OMP": "#ff9300", "Gg8": "#7a81ff"}

    selected_summary_fig = make_subplots(rows=3, cols=2, x_title="Dataset")

    genotypes = selected_avgs["Genotype"].unique()

    # pdb.set_trace()
    for genotype in genotypes:

        x_datasets = selected_avgs.loc[selected_avgs["Genotype"] == genotype][
            "Dataset"
        ]

        # mean trace peak amplitude
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Mean Trace Peak (pA)"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=1,
            col=1,
        )

        # Mean Onset Latency (ms)
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Mean Onset Latency (ms)"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=1,
            col=2,
        )

        # onset jitter
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Onset Jitter"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=2,
            col=1,
        )

        # mean trace onset latency
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Mean Trace Onset Latency (ms)"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=2,
            col=2,
        )

        # mean time to peak
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Mean Time to Peak (ms)"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=3,
            col=1,
        )

        # mean trace time to peak
        selected_summary_fig.add_trace(
            go.Box(
                x=x_datasets,
                y=selected_avgs.loc[selected_avgs["Genotype"] == genotype][
                    "Mean Trace Time to Peak (ms)"
                ].squeeze(),
                name=genotype,
                line=dict(color=genotype_color[genotype]),
                legendgroup=genotype,
            ),
            row=3,
            col=2,
        )

        # Update xaxis properties
        # selected_summary_fig.update_xaxes(autorange="reversed")
        # this defines the dataset order for x-axes
        dataset_order = [
            "non-injected",
            "3dpi",
            "5dpi",
            "dox_3dpi",
            "dox_4dpi",
            "dox_5dpi",
        ]

        selected_summary_fig.update_xaxes(
            categoryorder="array", categoryarray=dataset_order
        )

        # Update yaxis properties
        selected_summary_fig.update_yaxes(
            title_text="Mean Response Amplitude (pA)",
            row=1,
            col=1,
            autorange="reversed",
        )
        selected_summary_fig.update_yaxes(
            title_text="Mean Onset Latency (ms)", row=1, col=2
        )
        selected_summary_fig.update_yaxes(
            title_text="Mean Onset Jitter", row=2, col=1
        )
        selected_summary_fig.update_yaxes(
            title_text="Mean Trace Onset Latency (ms)", row=2, col=2
        )
        selected_summary_fig.update_yaxes(
            title_text="Mean Time to Peak (ms)", row=3, col=1
        )
        selected_summary_fig.update_yaxes(
            title_text="Mean Trace Time to Peak (ms)", row=3, col=2
        )

    selected_summary_fig.update_layout(
        # yaxis_title='Onset Latency (ms)',
        boxmode="group"  # group together boxes of the different traces for each value of x
    )

    # below is code from stack overflow to hide duplicate legends
    names = set()
    selected_summary_fig.for_each_trace(
        lambda trace: trace.update(showlegend=False)
        if (trace.name in names)
        else names.add(trace.name)
    )

    selected_summary_fig.update_layout(
        legend_title_text="Genotype",
        title_text=(
            "OMP vs. Gg8, {} ms onset latency threshold".format(threshold)
        ),
        title_x=0.5,
    )

    # selected_summary_fig.show()
    return selected_summary_fig


def save_selected_summary_fig(threshold, selected_summary_fig):

    html_filename = "{}_ms_threshold_datasets_summary.html".format(threshold)
    path = os.path.join(FileSettings.FIGURES_FOLDER, html_filename)

    selected_summary_fig.write_html(
        path, full_html=False, include_plotlyjs="cdn"
    )


def plot_response_counts(response_counts_dict):
    # response/no response is a trace
    thresholds = FileSettings.THRESHOLD_LIST.copy()
    dataset_order = [
        "non-injected/OMP",
        "non-injected/Gg8",
        "3dpi/Gg8",
        "5dpi/OMP",
        "5dpi/Gg8",
        "dox_3dpi/Gg8",
        "dox_4dpi/Gg8",
        "dox_5dpi/Gg8",
    ]

    response_colors = {"no response": "#A7BBC7", "response": "#293B5F"}
    response_counts_fig = make_subplots(
        rows=len(thresholds), cols=1, x_title="Dataset/Genotype"
    )

    # rearranges threshold for better plotting order
    thresholds.reverse()
    thresholds.insert(0, thresholds.pop(thresholds.index("nothresh")))
    for count, threshold in enumerate(thresholds):
        for response_type in response_counts_dict[threshold].keys():

            # pdb.set_trace()

            x_axis = response_counts_dict[threshold]["response"].columns
            response_counts_fig.add_trace(
                go.Bar(
                    x=x_axis,
                    y=response_counts_dict[threshold][response_type].squeeze(),
                    name=response_type,
                    marker_color=response_colors[response_type],
                    legendgroup=response_type,
                ),
                row=count + 1,
                col=1,
            )

            response_counts_fig.update_yaxes(
                title_text="{} ms latency cutoff cell count".format(threshold)
                if threshold != "nothresh"
                else "no onset latency cutoff cell count",
                row=count + 1,
                col=1,
            )

            response_counts_fig.update_layout(barmode="stack")

    # below is code from stack overflow to hide duplicate legends
    names = set()
    response_counts_fig.for_each_trace(
        lambda trace: trace.update(showlegend=False)
        if (trace.name in names)
        else names.add(trace.name)
    )

    response_counts_fig.update_xaxes(
        categoryorder="array", categoryarray=dataset_order
    )
    response_counts_fig.update_layout(
        legend_title_text="Cell Responses",
        title_text=(
            "Cell Responses by Onset Latency Cut-off".format(threshold)
        ),
        title_x=0.5,
    )
    response_counts_fig.show()
    return response_counts_fig


def save_response_counts_fig(response_counts_fig):

    html_filename = "all_cell_counts.html"
    path = os.path.join(FileSettings.FIGURES_FOLDER, html_filename)

    response_counts_fig.write_html(
        path, full_html=False, include_plotlyjs="cdn"
    )
