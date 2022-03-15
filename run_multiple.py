import os
import pandas as pd
from collections import defaultdict
from run_single_test import run_single

# from single_test import JaneCell
from aggregate_stats import GenotypeSummary
import pdb


def get_datasets(data_folder, ignored):
    dataset_list = [
        dataset
        for dataset in os.listdir(data_folder)
        if dataset not in ignored
    ]
    return dataset_list


def initialize_parameters(data_folder, ignored):
    # gets the list of datasets from file directory
    dataset_list = get_datasets(data_folder, ignored)
    # dataset_list = ["non-injected", "dox_5dpi"]

    # runs stats analysis for each dataset
    dataset_cell_counts = defaultdict(lambda: defaultdict(dict))
    threshold_list = [2, 4]

    return dataset_list, dataset_cell_counts, threshold_list


def run_dataset_analysis(data_folder, dataset, tables_folder):
    dataset_data_folder = os.path.join(data_folder, dataset)

    csvfile_name = "{}_sweep_info.csv".format(dataset)
    csvfile = os.path.join(tables_folder, dataset, csvfile_name)

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


def get_genotypes(tables_folder, dataset):
    # listing the genotypes in each dataset
    stats_folder = os.path.join(tables_folder, dataset)
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


def get_genotype_summary(
    dataset, genotypes_list, threshold_list, empty_count_dict
):
    for genotype in genotypes_list:
        genotype_summary = GenotypeSummary(
            dataset, genotype, tables_folder, figures_folder
        )
        genotype_summary.get_summary_stats()

        # counts_dict = {}  # one cell_counts df for each threshold
        for threshold in threshold_list:
            genotype_summary.set_latency_threshold(threshold)
            cell_counts = genotype_summary.count_cells()

            empty_count_dict[dataset][genotype][threshold] = cell_counts
            # genotype_summary.save_cell_counts()
            genotype_summary.calc_summary_avgs()
            genotype_summary.plot_averages()
            genotype_summary.save_summary_stats_fig()

    return empty_count_dict


def get_all_cell_counts(counts_dict, thresholds):
    """
    Takes the cell counts from all genotypes/datasets and compiles into one 
    df.
    """

    for threshold in thresholds:
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
    path = os.path.join(tables_folder, csv_filename)
    all_counts.to_csv(path, float_format="%8.4f")


def main():
    dataset_list, dataset_cell_counts, threshold_list = initialize_parameters(
        data_folder, ignored
    )

    for dataset_count, dataset in enumerate(dataset_list):
        print("***Starting analysis for {} dataset.***".format(dataset))
        # run_dataset_analysis(
        #     data_folder, dataset, tables_folder,
        # )
        genotypes_list = get_genotypes(tables_folder, dataset)
        dataset_cell_counts = get_genotype_summary(
            dataset, genotypes_list, threshold_list, dataset_cell_counts
        )

        print(
            "***Analysis for {} dataset done, #{}/{} datasets.***".format(
                dataset, dataset_count + 1, len(dataset_list)
            )
        )
    get_all_cell_counts(dataset_cell_counts, threshold_list)


if __name__ == "__main__":
    data_folder = "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data"
    ignored = {"esc_unusable"}
    tables_folder = (
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables"
    )
    figures_folder = (
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/figures"
    )

    main()

