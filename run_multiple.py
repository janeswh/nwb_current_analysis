import os
from run_single_test import run_single
from single_test import JaneCell
import pdb


def get_datasets(data_folder, ignored):
    dataset_list = [
        dataset
        for dataset in os.listdir(data_folder)
        if dataset not in ignored
    ]
    return dataset_list


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


def main():
    # gets the list of datasets from file directory
    dataset_list = get_datasets(data_folder, ignored)
    # dataset_list = ["dox_3dpi"]

    # runs stats analysis for each dataset
    for dataset_count, dataset in enumerate(dataset_list):
        print("***Starting analysis for {} dataset.***".format(dataset))
        run_dataset_analysis(
            data_folder, dataset, tables_folder,
        )
        print(
            "***Analysis for {} dataset done, #{}/{} datasets.***".format(
                dataset, dataset_count + 1, len(dataset_list)
            )
        )


if __name__ == "__main__":
    data_folder = "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data"
    ignored = {"esc_unusable"}
    tables_folder = (
        "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables"
    )

    main()

