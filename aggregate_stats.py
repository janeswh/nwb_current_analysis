import pandas as pd
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import sem
from collections import defaultdict
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

    def get_summary_avgs(self):
        # gets avg values for the selected stim condition
        self.selected_avgs = self.thresh_concat[
            (
                self.thresh_concat["Light Intensity"]
                == file_settings.selected_condition[0]
            )
            & (
                self.thresh_concat["Light Duration"]
                == file_settings.selected_condition[1]
            )
        ].copy()

        return self.selected_avgs

    def save_summary_avgs(self):
        csv_filename = "{}_{}_{}msthresh_summary_averages.csv".format(
            self.dataset, self.genotype, self.threshold
        )
        path = os.path.join(self.genotype_stats_folder, csv_filename)
        self.selected_avgs.to_csv(path, float_format="%8.4f", index=False)

    def count_monosyn_cells(self):
        # count the number of cells with monosynaptic response used for each
        # stimulus condition after thresholding with mean onset latency

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
        # pdb.set_trace()
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


def get_patched_counts(dataset_list):
    """
    Gets the number of files (i.e. # of cells patched) in all datasets
    using sweep_info csv files.
    """

    recorded_counts = defaultdict(lambda: defaultdict(dict))

    for dataset in dataset_list:
        csvfile_name = "{}_sweep_info.csv".format(dataset)
        csvfile = os.path.join(
            file_settings.tables_folder, dataset, csvfile_name
        )
        cell_list = pd.read_csv(csvfile, index_col=0)
        genotypes = cell_list["Genotype"].unique()

        for genotype in genotypes:
            genotype_count = len(
                cell_list.loc[cell_list["Genotype"] == genotype]
            )
            genotype_count = pd.DataFrame(
                {(dataset + "/" + genotype): genotype_count}, index=[0]
            )
            recorded_counts[dataset][genotype] = genotype_count

    return recorded_counts


def make_patched_counts_df(recorded_counts):
    all_patched = pd.DataFrame()

    for dataset in recorded_counts.keys():
        # use the genotypes found in dict
        for genotype in recorded_counts[dataset].keys():
            all_patched = pd.concat(
                [all_patched, recorded_counts[dataset][genotype]], axis=1,
            )
            all_patched.fillna(0, inplace=True)
            all_patched = all_patched.astype(int)

    return all_patched


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


def get_genotype_summary(dataset, genotypes_list, monosyn_count_dict):
    for genotype in genotypes_list:
        genotype_summary = GenotypeSummary(dataset, genotype)
        genotype_summary.get_summary_stats()

        # counts_dict = {}  # one cell_counts df for each threshold
        for threshold in file_settings.threshold_list:
            genotype_summary.set_latency_threshold(threshold)
            monosyn_cell_counts = genotype_summary.count_monosyn_cells()
            monosyn_count_dict[dataset][genotype][
                threshold
            ] = monosyn_cell_counts

            genotype_summary.get_summary_avgs()

            # genotype_summary.calc_summary_avgs()
            genotype_summary.save_summary_avgs()
            genotype_summary.plot_averages()
            genotype_summary.save_summary_stats_fig()

    return monosyn_count_dict


def collect_selected_averages(counts_dict):
    all_selected_averages = {}
    for threshold in file_settings.threshold_list:
        # pdb.set_trace()
        file_paths = []
        for dataset in counts_dict.keys():
            # use the genotypes found in dict
            for genotype in counts_dict[dataset].keys():
                for filename in os.listdir(
                    os.path.join(
                        file_settings.tables_folder, dataset, genotype
                    )
                ):
                    if filename.endswith(
                        "{}msthresh_summary_averages.csv".format(threshold)
                    ):
                        file_paths.append(
                            os.path.join(
                                file_settings.tables_folder,
                                dataset,
                                genotype,
                                filename,
                            )
                        )

        concat_avgs = pd.concat(
            [
                pd.read_csv(file, index_col=None, header=0)
                for file in file_paths
            ]
        )
        all_selected_averages[threshold] = concat_avgs

    return all_selected_averages


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
    path = os.path.join(file_settings.figures_folder, html_filename)

    selected_summary_fig.write_html(
        path, full_html=False, include_plotlyjs="cdn"
    )


def get_monosyn_cell_counts(monosyn_counts_dict):
    """
    Takes the monosynaptic cell counts from all genotypes/datasets and 
    compiles into one df.
    """

    for threshold in file_settings.threshold_list:
        all_counts = pd.DataFrame()
        # pdb.set_trace()

        for dataset in monosyn_counts_dict.keys():
            # use the genotypes found in dict
            for genotype in monosyn_counts_dict[dataset].keys():
                all_counts = pd.concat(
                    [
                        all_counts,
                        monosyn_counts_dict[dataset][genotype][threshold],
                    ],
                    axis=1,
                )
        all_counts.fillna(0, inplace=True)
        all_counts = all_counts.astype(int)

        save_cell_counts(all_counts, threshold)


def save_cell_counts(all_counts, threshold):
    """
    Takes the monosyn cell counts from specified threshold and saves as csv.
    """
    csv_filename = "{}ms_thresh_cell_counts.csv".format(threshold)
    path = os.path.join(file_settings.tables_folder, csv_filename)
    all_counts.to_csv(path, float_format="%8.4f")


def get_response_counts(all_counts, response_counts):
    """
    For each threshold, calculates:
    - % of cells that had a response / total # cells recorded in a dataset
    - % of cells that had a response / total # of cells that had a response 
    """
    # pdb.set_trace()
    all_counts.drop("dox_3dpi_ctrl/Gg8", axis=1, inplace=True)

    response_counts_dict = defaultdict(dict)

    for threshold in file_settings.threshold_list:

        csv_filename = "{}ms_thresh_cell_counts.csv".format(threshold)
        csvfile = os.path.join(file_settings.tables_folder, csv_filename)
        counts_df = pd.read_csv(csvfile, index_col=[0, 1])
        response_counts = pd.DataFrame(
            counts_df.loc[file_settings.selected_condition]
        ).T
        response_counts.index = [0]  # reset index for division
        no_response_counts = all_counts - response_counts

        response_counts_dict[threshold]["response"] = response_counts
        response_counts_dict[threshold]["no response"] = no_response_counts

    return response_counts_dict


def save_response_counts(all_patched, response_counts_dict):
    """
    Puts all the counts for total recorded cells and monosynaptic responding
    cells for each latency into one df, then export as csv.
    """
    all_counts = pd.DataFrame()

    all_patched = all_patched.copy()
    buffer_info = pd.DataFrame(
        {"Onset Latency Threshold (ms)": "N/A", "Count Type": "All patched",},
        index=[0],
    )
    all_patched = pd.concat([buffer_info, all_patched], axis=1)

    responses_dict = response_counts_dict.copy()
    for threshold in file_settings.threshold_list:
        for response in responses_dict[threshold].keys():
            info = pd.DataFrame(
                {
                    "Onset Latency Threshold (ms)": threshold,
                    "Count Type": response,
                },
                index=[0],
            )
            temp_df = pd.concat(
                [info, responses_dict[threshold][response]], axis=1
            )
            all_counts = pd.concat([all_counts, temp_df])

    all_counts = pd.concat([all_patched, all_counts])

    all_counts = all_counts[
        [
            "Onset Latency Threshold (ms)",
            "Count Type",
            "non-injected/OMP",
            "non-injected/Gg8",
            "3dpi/Gg8",
            "5dpi/OMP",
            "5dpi/Gg8",
            "dox_3dpi/Gg8",
            "dox_4dpi/Gg8",
            "dox_5dpi/Gg8",
        ]
    ]

    csv_filename = "monosynaptic_response_counts.csv".format(threshold)
    path = os.path.join(file_settings.tables_folder, csv_filename)
    all_counts.to_csv(path, float_format="%8.4f")

    return all_counts


def save_response_proportions(all_counts):

    response_proportions = pd.DataFrame()

    all_patched = all_counts[
        (all_counts["Count Type"] == "All patched")
        & (all_counts["Onset Latency Threshold (ms)"] == "N/A")
    ].copy()
    # drop str columns for division
    all_patched.drop(
        ["Onset Latency Threshold (ms)", "Count Type"], axis=1, inplace=True
    )

    nothresh_row = all_counts[
        (all_counts["Count Type"] == "response")
        & (all_counts["Onset Latency Threshold (ms)"] == "nothresh")
    ].copy()
    # drop str columns for division
    nothresh_row.drop(
        ["Onset Latency Threshold (ms)", "Count Type"], axis=1, inplace=True
    )

    for threshold in all_counts["Onset Latency Threshold (ms)"].unique():

        if not threshold in ["nothresh", "N/A"]:
            response_row = all_counts[
                (all_counts["Count Type"] == "response")
                & (all_counts["Onset Latency Threshold (ms)"] == threshold)
            ]
            #
            temp_response = response_row.drop(
                ["Onset Latency Threshold (ms)", "Count Type"], axis=1
            ).copy()

            response_allpatched = temp_response / all_patched
            response_allpatched[
                "Count Type"
            ] = "% responses out of all patched"
            response_allpatched["Onset Latency Threshold (ms)"] = threshold

            response_nothresh = temp_response / nothresh_row
            response_nothresh[
                "Count Type"
            ] = "% responses out of nothresh responses"
            response_nothresh["Onset Latency Threshold (ms)"] = threshold

            temp_proportions = pd.concat(
                [response_allpatched, response_nothresh]
            )

            response_proportions = pd.concat(
                [response_proportions, temp_proportions]
            )

    response_proportions = response_proportions[
        [
            "Onset Latency Threshold (ms)",
            "Count Type",
            "non-injected/OMP",
            "non-injected/Gg8",
            "3dpi/Gg8",
            "5dpi/OMP",
            "5dpi/Gg8",
            "dox_3dpi/Gg8",
            "dox_4dpi/Gg8",
            "dox_5dpi/Gg8",
        ]
    ]

    csv_filename = "monosynaptic_response_proportions.csv".format(threshold)
    path = os.path.join(file_settings.tables_folder, csv_filename)
    response_proportions.to_csv(path, float_format="%8.4f")


def plot_response_counts(response_counts_dict):
    # response/no response is a trace
    thresholds = file_settings.threshold_list.copy()
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
    path = os.path.join(file_settings.figures_folder, html_filename)

    response_counts_fig.write_html(
        path, full_html=False, include_plotlyjs="cdn"
    )


def do_cell_counts(monosyn_cell_counts, all_patched):
    """
    All the functions for counting cells with monosynaptic responses
    """
    # counts how many cells had monosynaptic responses after thresholding
    # with onset latency
    get_monosyn_cell_counts(monosyn_cell_counts)
    response_counts_dict = get_response_counts(
        all_patched, monosyn_cell_counts
    )

    # saves the counts and the calculated proportions for all datasets
    all_counts = save_response_counts(all_patched, response_counts_dict)
    save_response_proportions(all_counts)

    # plots and save the cell counts
    response_counts_fig = plot_response_counts(response_counts_dict)
    save_response_counts_fig(response_counts_fig)


def analyze_selected_condition(monosyn_cell_counts):
    """
    Pulls out the averages for each cell in the dataset for the selected
    stim condition (light intensity + duration).
    """
    all_selected_averages = collect_selected_averages(monosyn_cell_counts)

    for threshold in file_settings.threshold_list:
        selected_summary_fig = plot_selected_averages(
            threshold, all_selected_averages
        )
        save_selected_summary_fig(threshold, selected_summary_fig)

