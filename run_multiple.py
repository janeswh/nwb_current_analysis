import os
import pandas as pd
from collections import defaultdict
from run_single_test import run_single
from aggregate_stats import *
import file_settings

# from single_test import JaneCell
# from aggregate_stats import GenotypeSummary
import pdb


def get_datasets():
    dataset_list = [
        dataset
        for dataset in os.listdir(file_settings.data_folder)
        if dataset not in file_settings.ignored
    ]
    return dataset_list


def initialize_parameters():
    # gets the list of datasets from file directory
    dataset_list = get_datasets()
    # dataset_list = ["5dpi"]

    # runs stats analysis for each dataset
    dataset_cell_counts = defaultdict(lambda: defaultdict(dict))
    empty_selected_avgs = defaultdict(lambda: defaultdict(dict))

    return dataset_list, dataset_cell_counts, empty_selected_avgs


def run_dataset_analysis(dataset):
    dataset_data_folder = os.path.join(file_settings.data_folder, dataset)

    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(file_settings.tables_folder, dataset, csvfile_name)

    nwbfile_list = []

    for file in os.listdir(dataset_data_folder):
        if file.endswith(".nwb"):
            nwbfile_list.append(file)

    for file_count, nwbfile_name in enumerate(nwbfile_list):
        run_single(dataset, csvfile, nwbfile_name)
        print(
            "Analysis for {} done, #{}/{} cells".format(
                nwbfile_name, file_count + 1, len(nwbfile_list)
            )
        )


def main():
    (
        dataset_list,
        dataset_cell_counts,
        empty_selected_avgs,
    ) = initialize_parameters()

    for dataset_count, dataset in enumerate(dataset_list):
        print("***Starting analysis for {} dataset.***".format(dataset))
        # run_dataset_analysis(dataset)
        genotypes_list = get_genotypes(dataset)
        dataset_cell_counts, selected_avgs = get_genotype_summary(
            dataset, genotypes_list, dataset_cell_counts, empty_selected_avgs
        )

        print(
            "***Analysis for {} dataset done, #{}/{} datasets.***".format(
                dataset, dataset_count + 1, len(dataset_list)
            )
        )

    get_all_cell_counts(dataset_cell_counts)
    all_selected_averages = collect_selected_averages(dataset_cell_counts)

    for threshold in file_settings.threshold_list:
        selected_summary_fig = plot_selected_averages(
            threshold, all_selected_averages
        )
        save_selected_summary_fig(threshold, selected_summary_fig)


if __name__ == "__main__":

    main()

