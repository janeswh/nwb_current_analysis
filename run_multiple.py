import os
from run_single_test import run_single
from single_test import JaneCell

dataset = "non-injected"
data_folder = os.path.join(
    "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data", dataset
)

csvfile_name = "{}_sweep_info.csv".format(dataset)
csvfile = os.path.join(
    "/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/tables",
    dataset,
    csvfile_name,
)


nwbfile_list = []

for file in os.listdir(data_folder):
    if file.endswith(".nwb"):
        nwbfile_list.append(file)

for count, nwbfile_name in enumerate(nwbfile_list):
    run_single(dataset, csvfile, nwbfile_name)
    print(
        "Analysis for {} done, #{}/{} cells".format(
            nwbfile_name, count + 1, len(nwbfile_list)
        )
    )


# need to assign the same project dataset variable to all cells
# dataset variable determines path location to tables and figures


# phd_projects/MMZ_STC_dataset/tables
#   non-injected
#   5dpi
#   3 dpi, etc

