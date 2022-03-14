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


class GenotypeSummary(object):
    def __init__(self, dataset, genotype, tables_folder, figures_folder):
        self.dataset = dataset
        self.genotype = genotype
        self.genotype_stats_folder = os.path.join(
            tables_folder, dataset, genotype
        )
        self.genotype_figures_folder = os.path.join(
            figures_folder, dataset, genotype
        )
        self.concat_df = self.get_summary_stats()

    def get_summary_stats(self):
        # reads in all the summary stats values
        file_paths = []
        for filename in os.listdir(os.path.join(self.genotype_stats_folder)):
            if filename.endswith("response_stats.csv"):
                file_paths.append(
                    os.path.join(self.genotype_stats_folder, filename)
                )
        concat_df = pd.concat(
            [
                pd.read_csv(file, index_col=None, header=0)
                for file in file_paths
            ]
        )

        return concat_df

    def calc_summary_avgs(self):
        # this calculates stats grouped by stim conditions
        self.mean_df = self.concat_df.groupby(
            ["Light Intensity", "Light Duration"]
        ).mean()
        self.sem_df = self.concat_df.groupby(
            ["Light Intensity", "Light Duration"]
        ).sem()

        # this creates a list of unique stim condition combinations and allows iteration
        stim_conditions = (
            self.concat_df[["Light Intensity", "Light Duration"]]
            .drop_duplicates(ignore_index=True)
            .copy()
        )

        zip_conditions = zip(
            stim_conditions["Light Intensity"],
            stim_conditions["Light Duration"],
        )

        # this allows display of individual cell values matching each stim condition
        for intensity, duration in zip_conditions:
            indiv_cells_df = self.concat_df[
                (self.concat_df["Light Intensity"] == intensity)
                & (self.concat_df["Light Duration"] == duration)
            ].copy()

    def plot_averages(self):

        # plotting
        durations = self.concat_df["Light Duration"].unique()
        intensities = self.concat_df["Light Intensity"].unique()
        color_dict = {
            " 2 ms": "#7D1935",
            " 1 ms": "#B42B51",
            " 0.25 ms": "#E63E6D",
            " 0.01 ms": "#F892B9",
        }
        summary_stats_fig = make_subplots(
            rows=3, cols=2, x_title="Light Intensity (%)"
        )

        x_intensity = self.concat_df["Light Intensity"].unique()

        for count, duration in enumerate(durations):
            # mean trace peak amplitude
            summary_stats_fig.add_trace(
                go.Box(
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                    x=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.concat_df.loc[
                        self.concat_df["Light Duration"] == duration,
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
                self.dataset + " " + self.genotype + " summary values"
            ),
            title_x=0.5,
        )

        summary_stats_fig.show()
        self.summary_stats_fig = summary_stats_fig

    def save_summary_stats_fig(self):
        """
        Saves the summary stats plot as one html file
        """

        html_filename = "{}_summary_avgs.html".format(self.genotype)
        path = os.path.join(self.genotype_figures_folder, html_filename)

        self.summary_stats_fig.write_html(
            path, full_html=False, include_plotlyjs="cdn"
        )
