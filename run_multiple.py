import os

dataset = 'non-injected'
data_folder = os.path.join('/home/jhuang/Documents/phd_projects/MMZ_STC_dataset/data', dataset)
nwbfile_list = []

for file in os.listdir(data_folder):
    if file.endswith(".nwb"):
        nwbfile_list.append(file)
