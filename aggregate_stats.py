import pandas as pd
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import sem
import plotly.io as pio

pio.renderers.default = "browser"
import pdb
import glob

# def run_aggregate_stats(dataset):


dataset = "non-injected"
stats_folder = os.path.join(
    "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables", dataset
)

# loop to work on later

for genotype in os.listdir(stats_folder):
    file_paths = []
    for filename in os.listdir(os.path.join(stats_folder, genotype)):
        if filename.endswith("response_stats.csv"):
            file_paths.append(os.path.join(stats_folder, genotype, filename))
    concat_df = pd.concat(
        [pd.read_csv(file, index_col=None, header=0) for file in file_paths]
    )


# testing loop on one genotype

genotype = "Gg8"
test_stats_folder = os.path.join(
    "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
    dataset,
    genotype,
)

file_paths = []
for filename in os.listdir(os.path.join(test_stats_folder)):
    if filename.endswith("response_stats.csv"):
        file_paths.append(os.path.join(stats_folder, genotype, filename))
concat_df = pd.concat(
    [pd.read_csv(file, index_col=None, header=0) for file in file_paths]
)

# this calculates stats grouped by stim conditions
mean_df = concat_df.groupby(["Light Intensity", "Light Duration"]).mean()
sem_df = concat_df.groupby(["Light Intensity", "Light Duration"]).sem()

# this creates a list of unique stim condition combinations and allows iteration
stim_conditions = (
    concat_df[["Light Intensity", "Light Duration"]]
    .drop_duplicates(ignore_index=True)
    .copy()
)

zip_conditions = zip(
    stim_conditions["Light Intensity"], stim_conditions["Light Duration"]
)

# this allows display of individual cell values matching each stim condition
for intensity, duration in zip_conditions:
    indiv_cells_df = concat_df[
        (concat_df["Light Intensity"] == intensity)
        & (concat_df["Light Duration"] == duration)
    ].copy()


# plotting
durations = concat_df["Light Duration"].unique()
intensities = concat_df["Light Intensity"].unique()
color_dict = {
    " 2 ms": "#7D1935",
    " 1 ms": "#B42B51",
    " 0.25 ms": "#E63E6D",
    " 0.01 ms": "#F892B9",
}
curve_stats_fig = make_subplots(rows=3, cols=2, x_title="Light Intensity (%)")

x_intensity = concat_df["Light Intensity"].unique()

for count, duration in enumerate(durations):
    # mean trace peak amplitude
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration,
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
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration,
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
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Onset Jitter"],
            ].squeeze(),
            name=duration,
            line=dict(color=color_dict[duration]),
            legendgroup=duration,
        ),
        row=2,
        col=1,
    )

    # mean trace onset latency
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration,
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
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration,
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
    curve_stats_fig.add_trace(
        go.Box(
            x=concat_df.loc[
                concat_df["Light Duration"] == duration, ["Light Intensity"],
            ].squeeze(),
            y=concat_df.loc[
                concat_df["Light Duration"] == duration,
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
    # curve_stats_fig.update_xaxes(autorange="reversed")
    # this defines the intensities order for x-axes
    curve_stats_fig.update_xaxes(
        categoryorder="array", categoryarray=np.flip(intensities)
    )

    # Update yaxis properties
    curve_stats_fig.update_yaxes(
        title_text="Mean Response Amplitude (pA)",
        row=1,
        col=1,
        autorange="reversed",
    )
    curve_stats_fig.update_yaxes(
        title_text="Mean Onset Latency (ms)", row=1, col=2
    )
    curve_stats_fig.update_yaxes(title_text="Mean Onset Jitter", row=2, col=1)
    curve_stats_fig.update_yaxes(
        title_text="Mean Trace Onset Latency (ms)", row=2, col=2
    )
    curve_stats_fig.update_yaxes(
        title_text="Mean Time to Peak (ms)", row=3, col=1
    )
    curve_stats_fig.update_yaxes(
        title_text="Mean Trace Time to Peak (ms)", row=3, col=2
    )

curve_stats_fig.update_layout(
    # yaxis_title='Onset Latency (ms)',
    boxmode="group"  # group together boxes of the different traces for each value of x
)

# below is code from stack overflow to hide duplicate legends
names = set()
curve_stats_fig.for_each_trace(
    lambda trace: trace.update(showlegend=False)
    if (trace.name in names)
    else names.add(trace.name)
)

curve_stats_fig.update_layout(
    legend_title_text="Light Duration",
    title_text=(dataset + " " + genotype + " summary values"),
    title_x=0.5,
)


curve_stats_fig.show()

