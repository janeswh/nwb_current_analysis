import pandas as pd
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objects import Layout
from plotly.subplots import make_subplots
import plotly.io as pio

pio.kaleido.scope.default_scale = 5
pio.kaleido.scope.default_format = "png"
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


def plot_selected_averages(threshold, selected_avgs):

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
    thresholds[0] = "nothresh"
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


def plot_example_traces(
    trace1,
    trace2,
    type,
    same_genotype=None,
    condition_exception=False,
    small_yaxes=False,
):
    """
    Plots two mean traces on one plot for paper figures. Takes two traces as
    dfs and their type (e.g. ctrl vs. NBQX, Gg8 vs. OMP). Trace 1 is always
    control or OMP, Trace 2 is NBQX or Gg8. "genotype" indicates plots in which
    both cells are from the same genotype, for plotting color dictionary.
    """
    # sets background color to white
    layout = Layout(plot_bgcolor="rgba(0,0,0,0)")

    intensity, duration = FileSettings.SELECTED_CONDITION

    # sets Trace 1 values
    stim_columns = trace1.loc[:, ["Light Intensity", "Light Duration"]]
    trace1_to_plot = trace1.loc[
        :, 500.00:700.00
    ]  # only plots first 400-1000 ms

    trace1_to_plot_combined = pd.concat([stim_columns, trace1_to_plot], axis=1)

    trace1_y_toplot = trace1_to_plot_combined.loc[
        (trace1_to_plot_combined["Light Intensity"] == intensity)
        & (trace1_to_plot_combined["Light Duration"] == duration),
        500.00::,
    ].squeeze()

    # sets trace colors and Trace 2 values according to type of traces
    if type == "drug":
        if same_genotype == "OMP":
            type_names = ["OMP Control", "NBQX"]
            color = {"OMP Control": "#ff9300", "NBQX": "#EE251F"}
        elif same_genotype == "Gg8":
            type_names = ["Gg8 Control", "NBQX"]
            color = {"Gg8 Control": "#7a81ff", "NBQX": "#EE251F"}
        trace2_y_toplot = trace2.loc[
            500.00:700.00
        ].squeeze()  # only plots first 400-1000 ms

        trace2_x = trace2.loc[500.00:700.00].index

    elif type in {"genotypes", "OMP", "Gg8"}:
        if same_genotype is None:
            type_names = ["OMP", "Gg8"]
            color = {"OMP": "#ff9300", "Gg8": "#7a81ff"}
        elif same_genotype == "OMP":
            type_names = ["OMP cell 1", "OMP cell 2"]
            color = {"OMP cell 1": "#ff9300", "OMP cell 2": "#FBB85C"}
        elif same_genotype == "Gg8":
            type_names = ["Gg8 cell 1", "Gg8 cell 2"]
            color = {"Gg8 cell 1": "#7a81ff", "Gg8 cell 2": "#A4A8F9"}

        # chooses different condition for cell that doesn't have 100%, 1 ms
        if condition_exception is True:
            intensity = "50%"
            duration = " 1 ms"

        # color = {"OMP": "#ff9300", "Gg8": "#7a81ff"}
        trace2_to_plot = trace2.loc[
            :, 500.00:700.00
        ]  # only plots first 400-1000 ms

        trace2_to_plot_combined = pd.concat(
            [stim_columns, trace2_to_plot], axis=1
        )

        trace2_y_toplot = trace2_to_plot_combined.loc[
            (trace2_to_plot_combined["Light Intensity"] == intensity)
            & (trace2_to_plot_combined["Light Duration"] == duration),
            500.00::,
        ].squeeze()

        trace2_x = trace2_to_plot.columns

    example_traces_fig = go.Figure(layout=layout)
    example_traces_fig.add_trace(
        go.Scatter(
            x=trace1_to_plot.columns,
            y=trace1_y_toplot,
            name=type_names[0],
            mode="lines",
            line=dict(color=color[type_names[0]], width=4),
            # legendgroup=duration,
        ),
    )

    example_traces_fig.add_trace(
        go.Scatter(
            x=trace2_x,
            y=trace2_y_toplot,
            name=type_names[1],
            mode="lines",
            line=dict(color=color[type_names[1]], width=4),
            # legendgroup=duration,
        ),
    )

    # changes margin of plot to make room for light label
    # example_traces_fig.update_layout(margin=dict(t=150))
    # adds line for light stim
    example_traces_fig.add_shape(
        type="rect",
        x0=520,
        y0=5 if small_yaxes is True else 50,
        x1=521,
        y1=10 if small_yaxes is True else 100,
        line=dict(color="#33F7FF"),
        fillcolor="#33F7FF",
    )
    example_traces_fig.update_shapes(dict(xref="x", yref="y"))

    example_traces_fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        gridcolor="black",
        ticks="outside",
        tick0=520,
        dtick=10,
    )
    example_traces_fig.update_yaxes(
        showline=True, linewidth=1, gridcolor="black", linecolor="black",
    )

    if type == "genotypes":
        example_traces_fig.update_yaxes(range=[-700, 150])

    example_traces_fig.update_layout(
        font_family="Arial", legend=dict(font=dict(family="Arial", size=26))
    )

    example_traces_noaxes = go.Figure(example_traces_fig)
    example_traces_noaxes.update_xaxes(showgrid=False, visible=False)
    example_traces_noaxes.update_yaxes(showgrid=False, visible=False)

    return example_traces_fig, example_traces_noaxes


def save_example_traces_figs(
    axes, noaxes, type, same_genotype=None, cell_name=None
):
    """
    Saves the example traces figs as static png file
    """

    if not os.path.exists(FileSettings.PAPER_FIGURES_FOLDER):
        os.makedirs(FileSettings.PAPER_FIGURES_FOLDER)
    if type == "drug":
        axes_filename = "{}_NBQX_example_traces_axes.png".format(cell_name)
        noaxes_filename = "{}_NBQX_example_traces_noaxes.png".format(cell_name)
    elif type == "genotypes":
        if same_genotype is None:
            axes_filename = "OMPvGg8_traces_axes.png"
            noaxes_filename = "OMPvGg8_traces_noaxes.png"
        elif same_genotype == "OMP":
            axes_filename = "two_OMP_traces_axes.png"
            noaxes_filename = "two_OMP_traces_noaxes.png"
        elif same_genotype == "Gg8":
            axes_filename = "two_Gg8_traces_axes.png"
            noaxes_filename = "two_Gg8_traces_noaxes.png"

    axes.write_image(
        os.path.join(FileSettings.PAPER_FIGURES_FOLDER, axes_filename)
    )

    noaxes.write_image(
        os.path.join(FileSettings.PAPER_FIGURES_FOLDER, noaxes_filename)
    )


def plot_spike_sweep(trace):
    """
    Plots a single spike sweep to show STC physiology
    """
    layout = Layout(plot_bgcolor="rgba(0,0,0,0)")
    to_plot = trace[400:1600]
    spike_fig = go.Figure(layout=layout)
    spike_fig.add_trace(
        go.Scatter(
            x=trace.index,
            y=to_plot,
            # name=type_names[0],
            mode="lines",
            line=dict(color="#414145", width=2),
            # legendgroup=duration,
        ),
    )

    spike_fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        gridcolor="black",
        ticks="outside",
        tick0=400,
        dtick=100,
    )

    spike_fig.update_yaxes(
        showline=True, linewidth=1, gridcolor="black", linecolor="black",
    )

    spike_noaxes = go.Figure(spike_fig)
    spike_noaxes.update_xaxes(showgrid=False, visible=False)
    spike_noaxes.update_yaxes(showgrid=False, visible=False)

    pdb.set_trace()

    return spike_fig, spike_noaxes

