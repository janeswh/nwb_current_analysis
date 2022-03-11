import os
from run_single_test import run_single
from single_test import JaneCell


def main():
    ignored = {"esc_unusable"}
    data_folder = "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data"
    dataset_list = [
        dataset
        for dataset in os.listdir(data_folder)
        if dataset not in ignored
    ]

    print("***Starting analysis for {} dataset.***".format(dataset))

    for dataset_count, dataset in enumerate(dataset_list):
        dataset_data_folder = os.path.join(
            "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data", dataset
        )

        csvfile_name = "{}_sweep_info.csv".format(dataset)
        csvfile = os.path.join(
            "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
            dataset,
            csvfile_name,
        )

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

        print(
            "***Analysis for {} dataset done, #{}/{} datasets.***".format(
                dataset, dataset_count + 1, len(dataset_list)
            )
        )


if __name__ == "__main__":
    main()

