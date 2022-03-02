import os
import run_single_test
from single_test import JaneCell

dataset = 'non-injected'
data_folder = os.path.join('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data', dataset)
nwbfile_list = []

for file in os.listdir(data_folder):
    if file.endswith(".nwb"):
        nwbfile_list.append(file)


# need to assign the same project dataset variable to all cells 
# dataset variable determines path location to tables and figures


# phd_projects/MMZ_STC_dataset/tables
#   non-injected
#   5dpi
#   3 dpi, etc

