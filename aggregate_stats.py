import pandas as pd
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import sem
import plotly.io as pio
import file_settings

pio.renderers.default = "browser"
import pdb
import glob


class GenotypeSummary(object):
    def __init__(self, dataset, genotype):
        self.dataset = dataset
        self.genotype = genotype
        self.genotype_stats_folder = os.path.join(
            file_settings.tables_folder, dataset, genotype
        )
        self.genotype_figures_folder = os.path.join(
            file_settings.figures_folder, dataset, genotype
        )
        self.concat_df = self.get_summary_stats()

        self.threshold = None
        self.stim_conditions = None
        self.cell_counts = None

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

    def set_latency_threshold(self, threshold):
        # drops values from concat_df with mean onset latencies above threshold

        self.threshold = threshold
        if self.threshold == "nothresh":
            self.thresh_concat = self.concat_df.copy()
        else:
            self.thresh_concat = self.concat_df[
                self.concat_df["Mean Onset Latency (ms)"] < self.threshold
            ]

    def calc_summary_avgs(self):
        # this calculates stats grouped by stim conditions
        self.mean_df = self.thresh_concat.groupby(
            ["Light Intensity", "Light Duration"]
        ).mean()
        self.sem_df = self.thresh_concat.groupby(
            ["Light Intensity", "Light Duration"]
        ).sem()

    def save_summary_avgs(self):
        # saves avg values for the selected stim condition
        avgs_to_save = self.thresh_concat[
            (
                self.thresh_concat["Light Intensity"]
                == file_settings.selected_condition[0]
            )
            & (
                self.thresh_concat["Light Duration"]
                == file_settings.selected_condition[1]
            )
        ]

        # pdb.set_trace()

        csv_filename = "{}_{}_{}msthresh_summary_averages.csv".format(
            self.dataset, self.genotype, self.threshold
        )
        path = os.path.join(self.genotype_stats_folder, csv_filename)
        avgs_to_save.to_csv(path, float_format="%8.4f")

    def count_cells(self):
        # count the number of cells used for each stimulus condition after
        # thresholding with mean onset latency

        # this creates a list of unique stim condition combinations and allows iteration
        self.stim_conditions = (
            self.thresh_concat[["Light Intensity", "Light Duration"]]
            .drop_duplicates(ignore_index=True)
            .copy()
        )

        conditions_tuples = list(
            zip(
                self.stim_conditions["Light Intensity"],
                self.stim_conditions["Light Duration"],
            )
        )

        cell_counts = pd.DataFrame()
        counts_list = []
        # this allows display of individual cell values matching each stim condition
        for condition in conditions_tuples:
            indiv_cells_df = self.thresh_concat[
                (self.thresh_concat["Light Intensity"] == condition[0])
                & (self.thresh_concat["Light Duration"] == condition[1])
            ].copy()
            num_cells = len(indiv_cells_df)

            counts_list.append(num_cells)

        cell_counts = pd.DataFrame(
            {(self.dataset + "/" + self.genotype): counts_list}
        )
        cell_counts.index = pd.MultiIndex.from_tuples(conditions_tuples)

        return cell_counts

    def plot_averages(self):

        # plotting
        durations = self.thresh_concat["Light Duration"].unique()
        intensities = self.thresh_concat["Light Intensity"].unique()
        color_dict = {
            " 2 ms": "#7D1935",
            " 1 ms": "#B42B51",
            " 0.25 ms": "#E63E6D",
            " 0.01 ms": "#F892B9",
        }
        summary_stats_fig = make_subplots(
            rows=3, cols=2, x_title="Light Intensity (%)"
        )

        x_intensity = self.thresh_concat["Light Intensity"].unique()

        for count, duration in enumerate(durations):
            # mean trace peak amplitude
            summary_stats_fig.add_trace(
                go.Box(
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                    x=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
                        ["Light Intensity"],
                    ].squeeze(),
                    y=self.thresh_concat.loc[
                        self.thresh_concat["Light Duration"] == duration,
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
                self.dataset
                + " "
                + self.genotype
                + " summary values, "
                + str(self.threshold)
                + " mean onset latency threshold"
            ),
            title_x=0.5,
        )

        # summary_stats_fig.show()
        self.summary_stats_fig = summary_stats_fig

    def save_summary_stats_fig(self):

        html_filename = "{}_{}_threshold_summary_avgs.html".format(
            self.genotype, self.threshold
        )
        path = os.path.join(self.genotype_figures_folder, html_filename)

        self.summary_stats_fig.write_html(
            path, full_html=False, include_plotlyjs="cdn"
        )


def get_genotypes(dataset):
    # listing the genotypes in each dataset
    stats_folder = os.path.join(file_settings.tables_folder, dataset)
    genotypes_list = [
        genotype.name
        for genotype in os.scandir(stats_folder)
        if genotype.is_dir()
    ]

    if len(genotypes_list) == 0:
        print(
            "{} dataset has no cells with responses, no summary stats calculated.".format(
                dataset
            )
        )

    return genotypes_list


def get_genotype_summary(dataset, genotypes_list, empty_count_dict):
    for genotype in genotypes_list:
        genotype_summary = GenotypeSummary(dataset, genotype)
        genotype_summary.get_summary_stats()

        # counts_dict = {}  # one cell_counts df for each threshold
        for threshold in file_settings.threshold_list:
            genotype_summary.set_latency_threshold(threshold)
            cell_counts = genotype_summary.count_cells()

            empty_count_dict[dataset][genotype][threshold] = cell_counts
            # genotype_summary.save_cell_counts()
            genotype_summary.calc_summary_avgs()
            genotype_summary.save_summary_avgs()
            genotype_summary.plot_averages()
            genotype_summary.save_summary_stats_fig()

    return empty_count_dict


def get_all_cell_counts(counts_dict):
    """
    Takes the cell counts from all genotypes/datasets and compiles into one 
    df.
    """

    for threshold in file_settings.threshold_list:
        all_counts = pd.DataFrame()
        # pdb.set_trace()

        for dataset in counts_dict.keys():
            # use the genotypes found in dict
            for genotype in counts_dict[dataset].keys():
                all_counts = pd.concat(
                    [all_counts, counts_dict[dataset][genotype][threshold]],
                    axis=1,
                )
        all_counts.fillna(0, inplace=True)
        all_counts = all_counts.astype(int)
        save_cell_counts(all_counts, threshold)


def save_cell_counts(all_counts, threshold):
    """
    Takes the cell counts from specified threshold and saves as csv.
    """
    csv_filename = "{}ms_thresh_cell_counts.csv".format(threshold)
    path = os.path.join(file_settings.tables_folder, csv_filename)
    all_counts.to_csv(path, float_format="%8.4f")
